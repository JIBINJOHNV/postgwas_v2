import argparse
from argparse import Namespace
from pathlib import Path
import shutil
import sys
import polars as pl
import pytest
import yaml

from postgwas.config import load_configuration, load_module_configuration
from postgwas.config.models.modules.formatting import (
    FormattingStudyDesign,
    FormattingVcfFields,
)
from postgwas.core.variant_identifiers import inspect_bim_identifier_type
from postgwas.modules.formatting.cli import build_parser, main as formatter_main
from postgwas.modules.formatting.contracts import required_formats
from postgwas.modules.formatting.exporters.ldsc import export_ldsc
from postgwas.modules.formatting.exporters.mixer import export_mixer
from postgwas.modules.formatting.exporters.pred_ld import export_pred_ld
from postgwas.modules.formatting.service import run_formatter_direct
from postgwas.modules.formatting.table import (
    FormattingError,
    infer_study_design,
    select_variant_identifiers,
)


FIXTURE = Path(__file__).parent / "fixtures" / "formatting" / "harmonised.vcf"
TARGETS = ["magma", "gcta_gene", "susie", "finemap", "pred_ld", "ldsc", "mixer"]


def test_vcf_projection_uses_configured_canonical_names():
    projection = FormattingVcfFields.model_validate({"chromosome_value": "%CHROM"})
    assert projection.root == {"chromosome_value": "%CHROM"}


def test_pipeline_format_dependencies_come_from_formatter_configuration():
    config = load_configuration().modules.formatting.model_copy(
        update={"module_formats": {"custom_analysis": ["mixer"]}}
    )
    assert required_formats(config, ["custom_analysis"], ["ldsc"]) == ["ldsc", "mixer"]


def test_study_design_requirement_is_a_validated_formatter_contract():
    study_design = load_configuration().modules.formatting.study_design

    assert study_design.required_formats == ["ldsc", "mixer"]
    with pytest.raises(ValueError, match="ldsc and mixer exactly once"):
        FormattingStudyDesign.model_validate({
            "case_count_column": "N_CASE",
            "control_count_column": "N_CONTROL",
            "required_formats": ["ldsc"],
        })


def test_magma_susie_and_finemap_do_not_require_study_design_counts(
    tmp_path, monkeypatch, capsys,
):
    config = tmp_path / "formatting.yaml"
    config.write_text(
        "formats: [magma, susie, finemap]\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "CHROM": ["1"],
        "POS": [100],
        "ID": ["rs1"],
        "REF": ["A"],
        "ALT": ["G"],
        "BETA": [0.2],
        "SE": [0.1],
        "Z": [2.0],
        "LP": [8.0],
        "EAF": [0.2],
        "N": [1000.0],
        "NEFF": [900.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    output_directory = tmp_path / "results"
    results = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="QUANT",
        output_directory=str(output_directory),
        run_config=str(config),
        bcftools=sys.executable,
        overwrite=True,
    ))

    assert list(results) == ["magma", "susie", "finemap"]
    assert all(result["rows_out"] == 1 for result in results.values())
    assert pl.read_csv(
        results["magma"]["pval_file"], separator="\t",
    )["N_COL"].item() == 1000
    assert pl.read_csv(
        results["susie"]["susie_input"], separator="\t",
    )["NEF"].item() == 900
    assert pl.read_csv(
        results["finemap"]["finemap_input"], separator="\t",
    )["NEF"].item() == 900

    log_text = Path(results["magma"]["log_file"]).read_text(encoding="utf-8")
    assert "study_design" in log_text
    assert "reason=not_required_for_selected_formats" in log_text
    assert "Cannot infer the study type" not in log_text
    screen_text = capsys.readouterr().out
    assert "Study-design inference" in screen_text
    assert "not required for the selected formats" in screen_text
    resolved = yaml.safe_load(
        (output_directory / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert "study_design" not in resolved["modules"]["formatting"]


@pytest.mark.parametrize("target", ["ldsc", "mixer"])
def test_ldsc_and_mixer_still_require_study_design_counts(
    tmp_path, monkeypatch, target,
):
    config = tmp_path / (target + ".yaml")
    config.write_text("formats: [%s]\n" % target, encoding="utf-8")
    frame = pl.DataFrame({
        "CHROM": ["1"], "POS": [100], "ID": ["rs1"],
        "REF": ["A"], "ALT": ["G"], "Z": [2.0], "LP": [8.0],
        "EAF": [0.2], "INFO": [0.95], "NEFF": [900.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    with pytest.raises(
        FormattingError,
        match=r"canonical sample-count fields are missing: N_CASE, N_CONTROL",
    ):
        run_formatter_direct(Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(tmp_path / (target + "_results")),
            run_config=str(config),
            bcftools=sys.executable,
            overwrite=True,
        ))


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_one_vcf_extraction_creates_scientifically_consistent_outputs(tmp_path, capsys):
    args = Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=TARGETS,
        bcftools=shutil.which("bcftools"),
        overwrite=True,
    )

    result = run_formatter_direct(args)

    finemap = pl.read_csv(result["finemap"]["finemap_input"], separator="\t")
    assert finemap.columns == [
        "rsid", "chromosome", "position", "allele1", "allele2",
        "maf", "beta", "se", "NEF",
    ]
    rs1 = finemap.filter(pl.col("rsid") == "rs1").row(0, named=True)
    assert rs1["allele1"] == "G"  # ALT is the effect allele in GWAS-VCF
    assert rs1["allele2"] == "A"
    assert rs1["maf"] == pytest.approx(0.2)  # EAF=0.8 is not mislabeled as MAF
    assert finemap.height == 2  # missing rsID and SE=0 are independently invalid

    ldsc = pl.read_csv(result["ldsc"]["ldsc_file"], separator="\t")
    assert ldsc.columns == [
        "SNP", "A1", "A2", "Z", "P", "N_CAS", "N_CON", "FRQ", "INFO",
    ]
    assert ldsc.filter(pl.col("SNP") == "rs1")["A1"].item() == "G"
    assert ldsc.filter(pl.col("SNP") == "rs1")["N_CAS"].item() == 400
    assert ldsc.filter(pl.col("SNP") == "rs1")["N_CON"].item() == 600
    assert ldsc.filter(pl.col("SNP") == "rs3")["P"].item() == pytest.approx(1e-300)
    assert result["ldsc"]["sample_prev"] == pytest.approx(0.4)

    predld = pl.read_csv(
        Path(result["pred_ld"]["pred_ld_folder"]) / "STUDY_chr1_pred_ld_input.tsv",
        separator="\t",
    )
    assert predld.columns == [
        "snp", "chr", "pos", "A1", "A2", "beta", "SE",
        "NC", "SS", "AF", "LP", "SI",
    ]
    assert predld.filter(pl.col("snp") == "rs1")["NC"].item() == 600
    assert predld.filter(pl.col("snp") == "rs1")["SS"].item() == 1000

    susie = pl.read_csv(result["susie"]["susie_input"], separator="\t")
    assert susie.columns == ["SNP", "CHR", "BP", "REF", "ALT", "EZ", "LP", "NEF"]
    assert susie.select("REF", "ALT").row(0) == ("A", "G")
    magma = pl.read_csv(result["magma"]["pval_file"], separator="\t")
    assert magma.filter(pl.col("SNP") == "rs1")["N_COL"].item() == 1000
    magma_locations = pl.read_csv(
        result["magma"]["snp_loc_file"], separator="\t",
    )
    assert magma_locations.columns == ["SNP", "CHR", "BP", "REF", "ALT"]
    assert magma_locations.filter(pl.col("SNP") == "rs1").row(0) == (
        "rs1", 1, 100, "A", "G",
    )
    gcta = pl.read_csv(
        result["gcta_gene"]["summary_statistics_input_file"], separator="\t",
    )
    assert gcta.columns == ["SNP", "A1", "A2", "freq", "BETA", "SE", "P", "N"]
    assert gcta.filter(pl.col("SNP") == "rs1").row(0, named=True) == {
        "SNP": "rs1", "A1": "G", "A2": "A", "freq": 0.8,
        "BETA": 0.2, "SE": 0.1, "P": pytest.approx(1e-8), "N": 1000,
    }
    assert result["gcta_gene"]["sample_size_mode"] == "total_sample_size_from_FORMAT_SS"
    mixer = pl.read_csv(result["mixer"]["mixer_input"], separator="\t")
    assert mixer.columns == ["SNP", "CHR", "BP", "A1", "A2", "N", "Z"]
    assert mixer.filter(pl.col("SNP") == "rs1").row(0, named=True) == {
        "SNP": "rs1", "CHR": 1, "BP": 100, "A1": "G", "A2": "A",
        "N": 800.0, "Z": 2.0,
    }
    assert result["mixer"]["sample_size_mode"] == "effective_n_from_case_control_counts"
    assert (tmp_path / "logs" / "STUDY_formatter.log").is_file()
    resolved = yaml.safe_load(
        (tmp_path / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert list(resolved["modules"]) == ["formatting"]
    assert set(resolved["resources"]["executables"]) == {"bcftools"}
    assert resolved["run"]["dataset_id"] == "STUDY"
    assert resolved["run"]["output_directory"] == str(tmp_path)
    log_text = (tmp_path / "logs" / "STUDY_formatter.log").read_text(encoding="utf-8")
    assert "formatter_run status=COMPLETED" in log_text
    assert "DECIDE   study_design trait_type=binary" in log_text
    assert "metadata_header_used=false" in log_text
    screen_text = capsys.readouterr().out
    assert "Validated downstream input schemas" in screen_text
    assert "effective N" in screen_text
    assert "Inferred trait type" in screen_text
    assert "N_CASE present for 4/4 variants" in screen_text
    summary_lines = [
        line
        for line in screen_text.splitlines()
        if any(
            label in line
            for label in (
                "Dataset",
                "Input GWAS-VCF",
                "Inferred trait type",
                "magma",
                "Full log",
            )
        )
        and ":" in line
    ]
    assert len(summary_lines) == 5
    assert len({line.index(":") for line in summary_lines}) == 1


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_completed_formatter_outputs_resume_without_reextracting_vcf(
    tmp_path, monkeypatch,
):
    args = Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["gcta_gene"],
        bcftools=shutil.which("bcftools"),
        overwrite=True,
        resume=True,
    )
    first = run_formatter_direct(args)
    gcta = Path(first["gcta_gene"]["summary_statistics_input_file"])
    original_mtime = gcta.stat().st_mtime_ns
    assert (tmp_path / "run_metadata" / "formatter_completion.yaml").is_file()

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("resume re-read the GWAS-VCF"),
    )
    resumed = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["gcta_gene"],
        bcftools=shutil.which("bcftools"),
        resume=True,
    ))

    assert resumed["gcta_gene"]["resumed"] is True
    assert resumed["gcta_gene"]["summary_statistics_input_file"] == str(gcta)
    assert gcta.stat().st_mtime_ns == original_mtime


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_rebases_copied_run_to_current_output_directory(
    tmp_path, monkeypatch,
):
    original = tmp_path / "original"
    copied = tmp_path / "copied"
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        format=TARGETS,
        bcftools=shutil.which("bcftools"),
        resume=True,
    )
    first = run_formatter_direct(Namespace(
        **common, output_directory=str(original), overwrite=True,
    ))
    original_artifact = Path(
        first["gcta_gene"]["summary_statistics_input_file"]
    )
    shutil.copytree(original, copied)

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("copied resume re-read the GWAS-VCF"),
    )
    resumed = run_formatter_direct(Namespace(
        **common, output_directory=str(copied),
    ))

    current_artifact = (copied / "STUDY_gcta.ma").resolve()
    current_log = (copied / "logs" / "STUDY_formatter.log").resolve()
    result = resumed["gcta_gene"]
    assert Path(result["summary_statistics_input_file"]) == current_artifact
    assert Path(result["outputs"]["summary_statistics"]["path"]) == current_artifact
    assert Path(result["log_file"]) == current_log
    assert Path(result["summary_statistics_input_file"]) != original_artifact
    assert result["resumed"] is True
    assert Path(resumed["pred_ld"]["pred_ld_folder"]) == (
        copied / "pred_ld"
    ).resolve()
    returned_paths = {
        Path(value).resolve()
        for target_result in resumed.values()
        for value in _absolute_paths(target_result)
    }
    assert returned_paths
    assert all(
        path == copied.resolve() or copied.resolve() in path.parents
        for path in returned_paths
    )

    copied_manifest = yaml.safe_load(
        (copied / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert copied_manifest["results"]["gcta_gene"][
        "summary_statistics_input_file"
    ] == str(current_artifact)
    assert copied_manifest["results"]["gcta_gene"]["log_file"] == str(current_log)
    assert all(
        copied.resolve() in Path(record["path"]).parents
        for record in copied_manifest["outputs"]
    )


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_custom_output_resume_rebases_copied_artifact(tmp_path, monkeypatch):
    original = tmp_path / "original"
    copied = tmp_path / "copied"
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        bcftools=shutil.which("bcftools"),
        custom_output_file="study_custom.tsv",
        custom_columns={
            "id": str((tmp_path / "absolute-looking header").resolve()),
            "chr": "CHR",
        },
        resume=True,
    )
    first = run_formatter_direct(Namespace(
        **common, output_directory=str(original), overwrite=True,
    ))
    original_artifact = Path(first["custom"]["custom_output_file"])
    shutil.copytree(original, copied)

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("custom resume re-read the GWAS-VCF"),
    )
    resumed = run_formatter_direct(Namespace(
        **common, output_directory=str(copied),
    ))

    current_artifact = (copied / "study_custom.tsv").resolve()
    assert Path(resumed["custom"]["custom_output_file"]) == current_artifact
    assert current_artifact != original_artifact
    assert resumed["custom"]["resumed"] is True
    assert resumed["custom"]["field_roles"] == common["custom_columns"]
    manifest = yaml.safe_load(
        (copied / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["formats"] == []
    assert manifest["selected_targets"] == ["custom"]
    assert manifest["results"]["custom"]["custom_output_file"] == str(
        current_artifact
    )
    assert manifest["outputs"][0]["path"] == str(current_artifact)


def _absolute_paths(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from _absolute_paths(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _absolute_paths(child)
    elif isinstance(value, (str, Path)) and Path(value).is_absolute():
        yield value


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_rejects_manifest_filename_mismatch(tmp_path):
    original = tmp_path / "original"
    copied = tmp_path / "copied"
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        format=["gcta_gene"],
        bcftools=shutil.which("bcftools"),
        resume=True,
    )
    run_formatter_direct(Namespace(
        **common, output_directory=str(original), overwrite=True,
    ))
    shutil.copytree(original, copied)
    manifest_path = copied / "run_metadata" / "formatter_completion.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    wrong_artifact = str((original / "STUDY_wrong.ma").resolve())
    manifest["results"]["gcta_gene"][
        "summary_statistics_input_file"
    ] = wrong_artifact
    manifest["results"]["gcta_gene"]["outputs"]["summary_statistics"][
        "path"
    ] = wrong_artifact
    manifest_path.write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8",
    )

    with pytest.raises(
        FormattingError,
        match=r"artifact gcta_gene\.summary_statistics does not match.*STUDY_gcta\.ma",
    ):
        run_formatter_direct(Namespace(
            **common, output_directory=str(copied),
        ))

    log_text = (copied / "logs" / "STUDY_formatter.log").read_text(
        encoding="utf-8"
    )
    assert "FAILED" in log_text


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_does_not_fall_back_to_original_artifact(tmp_path):
    original = tmp_path / "original"
    copied = tmp_path / "copied"
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        format=["gcta_gene"],
        bcftools=shutil.which("bcftools"),
        resume=True,
    )
    first = run_formatter_direct(Namespace(
        **common, output_directory=str(original), overwrite=True,
    ))
    original_artifact = Path(
        first["gcta_gene"]["summary_statistics_input_file"]
    )
    shutil.copytree(original, copied)
    copied_artifact = copied / "STUDY_gcta.ma"
    copied_artifact.unlink()

    with pytest.raises(
        FormattingError,
        match=r"Resume file is missing or empty: .*copied/STUDY_gcta\.ma",
    ):
        run_formatter_direct(Namespace(
            **common, output_directory=str(copied),
        ))

    assert original_artifact.is_file()


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resolved_metadata_contains_only_selected_target(tmp_path):
    args = Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=["magma"], bcftools=shutil.which("bcftools"), overwrite=True,
    )

    run_formatter_direct(args)
    path = tmp_path / "run_metadata" / "resolved_config.yaml"
    resolved = yaml.safe_load(path.read_text(encoding="utf-8"))
    formatting = resolved["modules"]["formatting"]

    assert formatting["formats"] == ["magma"]
    assert list(formatting["exports"]) == ["magma"]
    assert "module_formats" not in formatting
    assert "format_order" not in formatting
    assert "chromosomes" not in formatting
    assert "study_design" not in formatting
    assert "mixer" not in formatting
    reloaded = load_configuration(path)
    assert reloaded.run.dataset_id == "STUDY"
    assert reloaded.run.output_directory == tmp_path
    assert reloaded.modules.formatting.formats == ["magma"]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_migrates_full_historical_configuration(tmp_path, monkeypatch):
    from postgwas.modules.formatting.resume import _configuration_values_digest

    args = Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=["magma"], bcftools=shutil.which("bcftools"), overwrite=True,
    )
    run_formatter_direct(args)

    resolved_path = tmp_path / "run_metadata" / "resolved_config.yaml"
    configuration = load_configuration(resolved_path)
    historical_module = configuration.modules.formatting.model_dump(mode="json")
    historical_module.pop("resolved_config", None)
    resolved = yaml.safe_load(resolved_path.read_text(encoding="utf-8"))
    resolved["modules"]["formatting"] = historical_module
    resolved_path.write_text(yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8")

    manifest_path = tmp_path / "run_metadata" / "formatter_completion.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("configuration_scope", None)
    manifest["configuration_sha256"] = _configuration_values_digest(
        historical_module,
        configuration.resources.executables.bcftools,
        ["magma"],
    )
    manifest_path.write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8",
    )

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("historical resume re-read the GWAS-VCF"),
    )
    resumed = run_formatter_direct(Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=["magma"], bcftools=shutil.which("bcftools"), resume=True,
    ))

    assert resumed["magma"]["resumed"] is True
    migrated = yaml.safe_load(resolved_path.read_text(encoding="utf-8"))
    assert list(migrated["modules"]["formatting"]["exports"]) == ["magma"]
    migrated_manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    assert migrated_manifest["configuration_scope"] == "selected_formatter_targets"
    assert migrated_manifest["configuration_sha256"] != manifest["configuration_sha256"]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_rejects_changed_input_vcf(tmp_path):
    vcf = tmp_path / "input.vcf"
    shutil.copyfile(FIXTURE, vcf)
    output = tmp_path / "output"
    common = dict(
        vcf=str(vcf), dataset_id="STUDY", output_directory=str(output),
        format=["gcta_gene"], bcftools=shutil.which("bcftools"), resume=True,
    )
    run_formatter_direct(Namespace(**common, overwrite=True))
    with vcf.open("a", encoding="utf-8") as handle:
        handle.write("\n")

    with pytest.raises(FormattingError, match="changed.*size mismatch"):
        run_formatter_direct(Namespace(**common))


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
@pytest.mark.parametrize("target", ["gcta_gene", "magma"])
def test_formatter_adopts_valid_completed_named_outputs_without_manifest(
    tmp_path, monkeypatch, target,
):
    common = dict(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=[target], bcftools=shutil.which("bcftools"),
    )
    run_formatter_direct(Namespace(**common, overwrite=True))
    manifest = tmp_path / "run_metadata" / "formatter_completion.yaml"
    manifest.unlink()

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("output adoption re-read the GWAS-VCF"),
    )
    resumed = run_formatter_direct(Namespace(**common))

    assert resumed[target]["adopted_existing_outputs"] is True
    assert manifest.is_file()
    resolved = yaml.safe_load(
        (tmp_path / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert "overwrite" not in resolved["modules"]["formatting"]
    assert resolved["run"]["overwrite"] is False


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_does_not_adopt_structurally_incomplete_magma_output(tmp_path):
    common = dict(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=["magma"], bcftools=shutil.which("bcftools"),
    )
    result = run_formatter_direct(Namespace(**common, overwrite=True))
    (tmp_path / "run_metadata" / "formatter_completion.yaml").unlink()
    with Path(result["magma"]["pval_file"]).open("a", encoding="utf-8") as handle:
        handle.write("broken\trecord\n")

    with pytest.raises(FormattingError, match="incomplete record"):
        run_formatter_direct(Namespace(**common))


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_overwrite_takes_precedence_over_default_resume(
    tmp_path, monkeypatch,
):
    args = Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY", output_directory=str(tmp_path),
        format=["magma"], bcftools=shutil.which("bcftools"), overwrite=True,
    )
    run_formatter_direct(args)

    from postgwas.modules.formatting import service

    monkeypatch.setattr(
        service,
        "resume_formatter_outputs",
        lambda **kwargs: pytest.fail("overwrite attempted to resume"),
    )
    rerun = run_formatter_direct(args)

    assert rerun["magma"].get("resumed") is not True


def test_formatter_rejects_colliding_outputs_before_vcf_extraction(
    tmp_path, monkeypatch,
):
    config = tmp_path / "formatting.yaml"
    config.write_text(
        "formats: [susie, finemap]\n"
        "exports:\n"
        "  susie:\n"
        "    output_file: '{dataset_id}_shared.tsv'\n"
        "  finemap:\n"
        "    output_file: '{dataset_id}_shared.tsv'\n",
        encoding="utf-8",
    )
    output = tmp_path / "results"
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("collision preflight read the VCF"),
    )

    with pytest.raises(FormattingError) as caught:
        run_formatter_direct(Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(output),
            run_config=str(config),
            bcftools=sys.executable,
            overwrite=True,
        ))

    message = str(caught.value)
    assert "Formatter output collision" in message
    assert "susie and finemap" in message
    assert str((output / "STUDY_shared.tsv").resolve()) in message
    assert not (output / "STUDY_shared.tsv").exists()
    assert not (output / "run_metadata" / "resolved_config.yaml").exists()
    log = output / "logs" / "STUDY_formatter.log"
    assert log.is_file()
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_formatter_rejects_colliding_partition_paths(tmp_path):
    from postgwas.modules.formatting.service import _validate_output_destinations

    module = load_configuration().modules.formatting
    pred_ld = module.exports["pred_ld"].model_copy(update={
        "partition_file": "{dataset_id}_pred_ld.tsv",
    })
    exports = dict(module.exports)
    exports["pred_ld"] = pred_ld
    module = module.model_copy(update={
        "chromosomes": ["1", "2"],
        "exports": exports,
    })

    with pytest.raises(
        FormattingError,
        match=(
            r"pred_ld\[partition=1\] and pred_ld\[partition=2\].*"
            r"Configure a unique"
        ),
    ):
        _validate_output_destinations(
            tmp_path,
            "STUDY",
            ["pred_ld"],
            module,
        )


@pytest.mark.parametrize(
    ("custom_output", "expected_label"),
    [
        ("STUDY_susie.tsv", "susie and custom"),
        ("logs/STUDY_formatter.log", "custom and formatter log"),
    ],
)
def test_custom_output_collision_fails_before_vcf_extraction(
    tmp_path, monkeypatch, custom_output, expected_label,
):
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("custom collision preflight read the VCF"),
    )

    with pytest.raises(FormattingError) as caught:
        run_formatter_direct(Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(tmp_path),
            format=["susie"],
            bcftools=sys.executable,
            custom_output_file=custom_output,
            custom_columns={"id": "SNP"},
            overwrite=True,
        ))

    assert "Formatter output collision" in str(caught.value)
    assert expected_label in str(caught.value)
    assert not (tmp_path / "STUDY_susie.tsv").exists()


def test_formatter_canonical_output_destinations_are_unique(tmp_path):
    from postgwas.modules.formatting.service import _validate_output_destinations

    module = load_configuration().modules.formatting

    destinations = _validate_output_destinations(
        tmp_path,
        "STUDY",
        TARGETS,
        module,
    )

    expected = sum(
        bool(module.exports[target].output_file)
        + len(module.exports[target].outputs)
        + (
            len(module.chromosomes)
            if module.exports[target].partition_file is not None
            else 0
        )
        for target in TARGETS
    )
    assert len(destinations) == expected
    assert len(set(destinations.values())) == expected


def test_pred_ld_reports_only_rows_written_for_configured_chromosomes(
    tmp_path, monkeypatch, capsys,
):
    config = tmp_path / "formatting.yaml"
    config.write_text(
        "formats: [pred_ld]\nchromosomes: ['1']\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "CHROM": ["1", "X"],
        "POS": [100, 200],
        "ID": ["rs1", "rs2"],
        "REF": ["A", "C"],
        "ALT": ["G", "T"],
        "BETA": [0.2, -0.1],
        "SE": [0.1, 0.2],
        "EAF": [0.2, 0.3],
        "INFO": [0.95, 0.90],
        "LP": [8.0, 6.0],
        "N": [1000.0, 1000.0],
        "N_CASE": [400.0, 400.0],
        "N_CONTROL": [600.0, 600.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "results"),
        run_config=str(config),
        bcftools=sys.executable,
        overwrite=True,
    ))["pred_ld"]

    assert result["rows_in"] == 2
    assert result["rows_out"] == 1
    assert result["rows_excluded_unconfigured_chromosome"] == 1
    assert result["excluded_chromosomes"] == {"X": 1}
    assert result["schema_rows_excluded"] == 0
    assert result["rows_excluded"] == 1
    output = pl.read_csv(result["files"][0], separator="\t")
    assert output["snp"].to_list() == ["rs1"]

    log_text = Path(result["log_file"]).read_text(encoding="utf-8")
    assert "excluded_unconfigured_chromosomes" in log_text
    assert "chromosomes='{X -> 1}'" in log_text
    assert "reason=unconfigured_chromosome" in log_text
    screen_text = capsys.readouterr().out
    assert "1 outside configured" in screen_text
    assert "chromosomes" in screen_text


def test_pred_ld_fails_before_writing_when_no_configured_chromosomes_remain(
    tmp_path,
):
    frame = pl.DataFrame({
        "SNP": ["rsX"], "CHROM": ["X"], "POS": [100],
        "ALT": ["G"], "REF": ["A"], "BETA": [0.2], "SE": [0.1],
        "N_CONTROL": [600.0], "N": [1000.0], "EAF": [0.2],
        "LP": [8.0], "INFO": [0.95],
    })
    config = load_configuration().modules.formatting.model_copy(
        update={"chromosomes": ["1"]}
    )

    with pytest.raises(FormattingError) as caught:
        export_pred_ld(frame, tmp_path, "STUDY", config, overwrite=False)

    message = str(caught.value)
    assert "No PRED-LD variants remain on the configured chromosomes: 1" in message
    assert "Excluded unconfigured chromosome counts: X=1" in message
    assert not (tmp_path / "pred_ld").exists()


def test_formatter_help_explains_formats_and_configuration():
    help_text = build_parser().format_help()
    assert "The VCF is read once" in help_text
    assert "pred_ld" in help_text
    assert "mixer" in help_text
    assert "modules.formatting.formats" in help_text
    assert "--variant-id-type {rsid,unique}" in help_text
    assert "Per-target YAML settings" in help_text
    assert "--run-config PATH" in help_text
    assert "--custom-output FILE" in help_text
    assert "--id NAME" in help_text
    assert "--n-case NAME" in help_text
    assert "without editing YAML" in help_text


def test_bare_formatter_command_displays_the_same_help_as_help_flag(capsys):
    assert formatter_main([]) == 0
    bare_help = capsys.readouterr().out

    with pytest.raises(SystemExit) as caught:
        formatter_main(["--help"])
    explicit_help = capsys.readouterr().out

    assert caught.value.code == 0
    assert bare_help == explicit_help
    assert "custom output" in bare_help.lower()


def test_formatter_options_without_vcf_still_fail_validation():
    with pytest.raises(SystemExit) as caught:
        formatter_main(["--format", "magma"])

    assert caught.value.code == 2


def test_formatter_configurable_cli_actions_do_not_own_defaults():
    parser = build_parser()
    for destination in (
        "format", "run_config", "resume", "overwrite", "bcftools", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed", "variant_id_type",
        "custom_output_file", "custom_columns",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS


def test_custom_cli_preserves_requested_column_order():
    args = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--custom-output", "study.tsv",
        "--p", "PVALUE",
        "--id", "MARKER",
        "--alt", "EFFECT",
        "--ref", "OTHER",
        "--lp", "LOGP",
        "--maf", "MAF",
    ])

    assert args.custom_columns == {
        "p": "PVALUE",
        "id": "MARKER",
        "alt": "EFFECT",
        "ref": "OTHER",
        "lp": "LOGP",
        "maf": "MAF",
    }


@pytest.mark.parametrize(
    ("attributes", "message"),
    [
        (
            {"custom_output_file": "study.tsv", "custom_columns": {"chr": "CHR"}},
            r"custom output requires --id NAME",
        ),
        (
            {"custom_columns": {"id": "SNP"}},
            r"custom output requires --custom-output",
        ),
        (
            {
                "custom_output_file": "study.tsv",
                "custom_columns": {"id": "SAME", "chr": "SAME"},
            },
            r"column names must be unique",
        ),
    ],
)
def test_custom_cli_rejects_incomplete_or_ambiguous_schema(attributes, message):
    from postgwas.modules.formatting.service import _resolved_configuration

    with pytest.raises(Exception, match=message):
        _resolved_configuration(Namespace(**attributes))


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_custom_output_is_additive_and_uses_requested_scientific_semantics(tmp_path):
    args = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path),
        "--bcftools", shutil.which("bcftools"),
        "--format", "magma",
        "--custom-output", "study_custom.tsv",
        "--p", "PVALUE",
        "--id", "MARKER",
        "--maf", "MAF",
        "--eaf", "EAF",
        "--alt", "EFFECT",
        "--ref", "OTHER",
        "--lp", "LOGP",
        "--variant-id-type", "unique",
        "--overwrite",
    ])

    result = run_formatter_direct(args)

    assert list(result) == ["magma", "custom"]
    assert Path(result["magma"]["pval_file"]).is_file()
    custom = pl.read_csv(result["custom"]["custom_output_file"], separator="\t")
    assert custom.columns == [
        "PVALUE", "MARKER", "MAF", "EAF", "EFFECT", "OTHER", "LOGP",
    ]
    assert custom.height == 4
    rs1 = custom.filter(pl.col("MARKER") == "1_100_A_G").row(0, named=True)
    assert rs1 == {
        "PVALUE": pytest.approx(1e-8),
        "MARKER": "1_100_A_G",
        "MAF": pytest.approx(0.2),
        "EAF": pytest.approx(0.8),
        "EFFECT": "G",
        "OTHER": "A",
        "LOGP": pytest.approx(8.0),
    }
    bounded = custom.filter(pl.col("MARKER") == "2_300_G_A")["PVALUE"].item()
    assert bounded == pytest.approx(1e-300)
    assert result["custom"]["p_values_bounded"] == 1
    assert result["custom"]["variant_id_type"] == "unique"
    log_text = Path(result["custom"]["log_file"]).read_text(encoding="utf-8")
    assert "target=custom" in log_text
    assert "study_design reason=not_required_for_selected_formats" in log_text

    resolved = yaml.safe_load(
        (tmp_path / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )["modules"]["formatting"]
    assert list(resolved["exports"]) == ["magma"]
    assert resolved["custom_output"]["columns"] == args.custom_columns


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_custom_id_supports_rsid_and_counts_unusable_identifiers(tmp_path):
    args = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--dataset-id", "STUDY",
        "--output-directory", str(tmp_path),
        "--bcftools", shutil.which("bcftools"),
        "--custom-output", "rsids.tsv",
        "--id", "src",
        "--chr", "CHR",
    ])

    result = run_formatter_direct(args)["custom"]
    output = pl.read_csv(result["custom_output_file"], separator="\t")

    assert output.columns == ["src", "CHR"]
    assert output["src"].to_list() == ["rs1", "rs3", "rs4"]
    assert result["rows_in"] == 4
    assert result["rows_out"] == 3
    assert result["identifier_rows_excluded"] == 1
    assert result["rows_excluded"] == 1
    manifest = yaml.safe_load(
        (tmp_path / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert [record["path"] for record in manifest["outputs"]] == [
        result["custom_output_file"]
    ]


def test_custom_output_excludes_and_reports_missing_requested_values(
    tmp_path, monkeypatch,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "1", "1"],
        "POS": [100, 200, 300, 400],
        "ID": ["rs1", "rs2", "rs3", "rs4"],
        "REF": ["A", "C", "G", "T"],
        "ALT": ["G", "T", "A", "C"],
        "BETA": [0.2, None, 0.3, 0.4],
        "SE": [0.1, 0.2, 0.0, 0.3],
        "EAF": [0.8, 0.2, 0.3, 1.2],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        bcftools=sys.executable,
        custom_output_file="study.tsv",
        custom_columns={
            "id": "SNP", "beta": "BETA", "se": "SE", "maf": "MAF",
        },
        overwrite=True,
    ))["custom"]

    output = pl.read_csv(result["custom_output_file"], separator="\t")
    assert output.to_dicts() == [{
        "SNP": "rs1", "BETA": 0.2, "SE": 0.1, "MAF": pytest.approx(0.2),
    }]
    assert result["rows_in"] == 4
    assert result["rows_out"] == 1
    assert result["schema_rows_excluded"] == 3
    assert result["rows_excluded"] == 3
    assert "missing_or_invalid_required_value" in Path(
        result["log_file"]
    ).read_text(encoding="utf-8")


def test_identifier_selection_is_general_and_configuration_driven():
    frame = pl.DataFrame({
        "CHROM": ["1", "2"],
        "POS": [100, 200],
        "ID": ["alias;rs10", "."],
        "REF": ["A", "C"],
        "ALT": ["G", "T"],
    })
    module = load_configuration().modules.formatting

    rsids, rsid_qc = select_variant_identifiers(frame, module, "rsid")
    unique, unique_qc = select_variant_identifiers(frame, module, "unique")

    assert rsids["SNP"].to_list() == ["rs10"]
    assert rsid_qc["identifier_rows_excluded"] == 1
    assert unique["SNP"].to_list() == ["1_100_A_G", "2_200_C_T"]
    assert unique_qc["identifier_rows_excluded"] == 0


def test_identifier_configuration_rejects_invalid_regex_and_template(tmp_path):
    invalid_regex = tmp_path / "invalid_regex.yaml"
    invalid_regex.write_text(
        "variant_identifiers:\n  rsid_extraction_pattern: '(rs[0-9]+)(extra)'\n",
        encoding="utf-8",
    )
    with pytest.raises(Exception, match="exactly one capturing group"):
        load_module_configuration("formatting", invalid_regex)

    invalid_template = tmp_path / "invalid_template.yaml"
    invalid_template.write_text(
        "variant_identifiers:\n"
        "  unique_id_template: '{chromosome}_{position}'\n",
        encoding="utf-8",
    )
    with pytest.raises(
        Exception,
        match="chromosome, position, reference_allele, and alternate_allele",
    ):
        load_module_configuration("formatting", invalid_template)


def test_cli_identifier_type_overrides_per_target_yaml(tmp_path):
    from postgwas.modules.formatting.service import _resolved_configuration

    config = tmp_path / "formatting.yaml"
    config.write_text(
        "variant_identifiers:\n"
        "  target_types: {magma: unique, ldsc: unique}\n",
        encoding="utf-8",
    )

    resolved = _resolved_configuration(Namespace(
        run_config=str(config),
        variant_id_type="rsid",
    )).modules.formatting.variant_identifiers

    assert resolved.default_type == "rsid"
    assert resolved.target_types == {}


@pytest.mark.parametrize(
    ("rows", "expected"),
    [
        (["1 rs1 0 100 G A", "1 rs2 0 200 T C"], "rsid"),
        (["chr1 1_100_A_G 0 100 G A", "1 1_200_C_T 0 200 T C"], "unique"),
    ],
)
def test_complete_bim_scan_detects_supported_identifier_type(tmp_path, rows, expected):
    bim = tmp_path / "reference.bim"
    bim.write_text("\n".join(rows) + "\n", encoding="utf-8")
    configuration = load_configuration()
    formatting = configuration.modules.formatting.variant_identifiers
    magma = configuration.modules.magma.input

    summary = inspect_bim_identifier_type(
        bim,
        column_roles=list(magma.bim_columns),
        delimiter_pattern=magma.table_delimiter_pattern,
        rsid_pattern=formatting.rsid_pattern,
        unique_id_template=formatting.unique_id_template,
        chromosome_prefix_pattern=(
            configuration.modules.formatting.chromosome_labels.prefix_pattern
        ),
        chromosome_aliases=configuration.modules.formatting.chromosome_labels.aliases,
    )

    assert summary.identifier_type == expected
    assert summary.variants == 2


def test_bim_scan_rejects_mixed_identifier_types(tmp_path):
    bim = tmp_path / "mixed.bim"
    bim.write_text(
        "1 rs1 0 100 G A\n1 1_200_C_T 0 200 T C\n",
        encoding="utf-8",
    )
    configuration = load_configuration()
    formatting = configuration.modules.formatting.variant_identifiers
    magma = configuration.modules.magma.input

    with pytest.raises(ValueError, match="mixed or unsupported"):
        inspect_bim_identifier_type(
            bim,
            column_roles=list(magma.bim_columns),
            delimiter_pattern=magma.table_delimiter_pattern,
            rsid_pattern=formatting.rsid_pattern,
            unique_id_template=formatting.unique_id_template,
            chromosome_prefix_pattern=(
                configuration.modules.formatting.chromosome_labels.prefix_pattern
            ),
            chromosome_aliases=(
                configuration.modules.formatting.chromosome_labels.aliases
            ),
        )


@pytest.mark.parametrize("analysis_module", ("magma", "pops", "magmacovar", "flames"))
def test_magma_analysis_paths_infer_only_the_magma_formatter_identifier_type(
    tmp_path, monkeypatch, analysis_module,
):
    from postgwas.pipeline.runners import run_formatter_runner

    reference = tmp_path / "reference"
    reference.with_suffix(".bim").write_text(
        "1 1_100_A_G 0 100 G A\n1 1_200_C_T 0 200 T C\n",
        encoding="utf-8",
    )
    captured = {}

    def fake_formatter(args, ctx):
        captured["target_types"] = dict(args.variant_id_types)
        captured["observation"] = dict(args.variant_id_observations["magma"])
        captured["formats"] = list(args.format)
        return {"magma": {"ok": True}}

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        fake_formatter,
    )
    output = tmp_path / "output"
    args = Namespace(
        modules=[analysis_module],
        vcf=str(FIXTURE),
        magma_ld_reference=str(reference),
        output_directory=str(output),
        dataset_id="STUDY",
        bcftools=shutil.which("bcftools") or sys.executable,
        _step_num="01",
    )

    result = run_formatter_runner(args, {})

    assert result == {"magma": {"ok": True}}
    assert captured["target_types"] == {"magma": "unique"}
    assert "magma" in captured["formats"]
    assert captured["observation"]["unique_ids"] == 2
    assert args.output_directory == str(output)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_per_target_identifier_types_do_not_change_other_exports(tmp_path):
    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["magma", "ldsc"],
        variant_id_types={"magma": "unique", "ldsc": "rsid", "mixer": "unique"},
        bcftools=shutil.which("bcftools"),
        overwrite=True,
    ))

    magma = pl.read_csv(result["magma"]["pval_file"], separator="\t")
    ldsc = pl.read_csv(result["ldsc"]["ldsc_file"], separator="\t")
    assert magma["SNP"].to_list() == [
        "1_100_A_G", "1_200_C_T", "2_300_G_A", "2_400_T_C",
    ]
    assert ldsc["SNP"].to_list() == ["rs1", "rs3", "rs4"]
    assert result["magma"]["variant_id_type"] == "unique"
    assert result["ldsc"]["variant_id_type"] == "rsid"
    resolved = yaml.safe_load(
        (tmp_path / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert resolved["modules"]["formatting"]["variant_identifiers"][
        "target_types"
    ] == {"magma": "unique", "ldsc": "rsid"}


def test_exported_module_only_configuration_is_accepted(tmp_path):
    config = tmp_path / "formatting.yaml"
    config.write_text(
        "enabled: true\nformats: [magma]\n"
        "chromosomes: ['1']\nminimum_p_value: 1.0e-250\n",
        encoding="utf-8",
    )
    args = Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(config), bcftools=shutil.which("bcftools") or "bcftools",
    )
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools is required")
    result = run_formatter_direct(args)
    assert list(result) == ["magma"]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_ldsc_column_names_and_output_path_are_configuration_driven(tmp_path):
    config = tmp_path / "formatting.yaml"
    config.write_text(
        "enabled: true\nformats: [ldsc]\n"
        "exports:\n"
        "  ldsc:\n"
        "    output_file: '{dataset_id}_custom_ldsc.tsv'\n"
        "    columns: {SNP: MARKER, ALT: EFFECT, REF: OTHER, Z: STAT, LP: PVALUE}\n"
        "    trait_columns:\n"
        "      binary: {N_CASE: CASES, N_CONTROL: CONTROLS}\n"
        "      quantitative: {N_CONTROL: TOTAL}\n"
        "    trailing_columns: {EAF: FREQUENCY, INFO: QUALITY}\n",
        encoding="utf-8",
    )
    args = Namespace(
        vcf=str(FIXTURE), dataset_id="STUDY",
        output_directory=str(tmp_path / "out"), run_config=str(config),
        bcftools=shutil.which("bcftools"),
    )

    result = run_formatter_direct(args)
    output = Path(result["ldsc"]["ldsc_file"])
    assert output.name == "STUDY_custom_ldsc.tsv"
    assert pl.read_csv(output, separator="\t").columns == [
        "MARKER", "EFFECT", "OTHER", "STAT", "PVALUE",
        "CASES", "CONTROLS", "FREQUENCY", "QUALITY",
    ]


def test_failure_still_writes_the_canonical_log(tmp_path):
    args = Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["magma"],
        bcftools="bcftools-that-does-not-exist",
    )
    with pytest.raises(RuntimeError, match="bcftools executable was not found"):
        run_formatter_direct(args)
    log = tmp_path / "logs" / "STUDY_formatter.log"
    assert log.is_file()
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_configuration_failure_still_writes_the_canonical_log(tmp_path):
    config = tmp_path / "invalid_formatting.yaml"
    config.write_text("formats: [unsupported_tool]\n", encoding="utf-8")
    args = Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(config),
    )

    with pytest.raises(Exception, match="modules.formatting.formats"):
        run_formatter_direct(args)

    log = tmp_path / "out" / "logs" / "STUDY_formatter.log"
    assert log.is_file()
    assert "Formatter configuration failed" in log.read_text(encoding="utf-8")


def test_ldsc_quantitative_schema_uses_total_sample_size(tmp_path):
    frame = pl.DataFrame({
        "SNP": ["rs1"], "ALT": ["G"], "REF": ["A"],
        "Z": [2.0], "LP": [8.0], "N": [9999.0],
        "N_CASE": [None], "N_CONTROL": [1000.0],
        "EAF": [0.2], "INFO": [0.95],
    })
    config = load_configuration().modules.formatting

    result = export_ldsc(frame, tmp_path, "QUANT", config, overwrite=False)
    output = pl.read_csv(result["ldsc_file"], separator="\t")

    assert output.columns == ["SNP", "A1", "A2", "Z", "P", "N", "FRQ", "INFO"]
    assert output["N"].item() == 1000
    assert result["trait_type"] == "quantitative"
    assert result["sample_size_mode"] == "quantitative_n_from_nco"
    assert result["sample_prev"] is None


def test_formatter_table_serialisation_uses_runtime_configuration(tmp_path):
    frame = pl.DataFrame({
        "SNP": ["rs1"], "ALT": ["G"], "REF": ["A"],
        "Z": [2.0], "LP": [8.0], "N": [1000.0],
        "N_CASE": [None], "N_CONTROL": [1000.0],
        "EAF": [0.2], "INFO": [None],
    })
    config = load_configuration().modules.formatting
    runtime = config.runtime.model_copy(update={
        "table_delimiter": ",",
        "output_null_value": "MISSING",
        "atomic_output_suffix": ".partial",
        "atomic_uncompressed_suffix": ".plain",
    })
    config = config.model_copy(update={"runtime": runtime})

    result = export_ldsc(frame, tmp_path, "CUSTOM", config, overwrite=False)
    text = Path(result["ldsc_file"]).read_text(encoding="utf-8")

    assert text.splitlines()[0] == "SNP,A1,A2,Z,P,N,FRQ,INFO"
    assert text.splitlines()[1].endswith(",MISSING")
    assert not list(tmp_path.glob("*.partial"))
    assert not list(tmp_path.glob("*.plain"))


def test_study_design_uses_values_not_vcf_metadata():
    quantitative = pl.DataFrame({
        "N_CASE": [None, None],
        "N_CONTROL": [1000.0, 950.0],
    })
    binary_with_one_missing_case = pl.DataFrame({
        "N_CASE": [400.0, None],
        "N_CONTROL": [600.0, 550.0],
    })

    assert infer_study_design(
        quantitative, "N_CASE", "N_CONTROL",
    ).trait_type == "quantitative"
    assert infer_study_design(
        binary_with_one_missing_case, "N_CASE", "N_CONTROL",
    ).trait_type == "binary"


def test_study_design_requires_nco_values():
    frame = pl.DataFrame({
        "N_CASE": [None, None],
        "N_CONTROL": [None, None],
    })
    with pytest.raises(FormattingError, match="N_CONTROL.*missing for every variant"):
        infer_study_design(frame, "N_CASE", "N_CONTROL")


def test_mixer_quantitative_n_uses_harmonised_neff(tmp_path):
    frame = pl.DataFrame({
        "SNP": ["rs1", "rs2"], "CHROM": ["1", "1"], "POS": [1, 2],
        "ALT": ["A", "C"], "REF": ["G", "T"], "Z": [1.0, -1.0],
        "NEFF": [1000.0, 300.0], "INFO": [0.9, 0.9],
        "N_CASE": [None, None], "N_CONTROL": [1000.0, 300.0],
    })
    config = load_configuration().modules.formatting

    result = export_mixer(frame, tmp_path, "QUANT", config, overwrite=False)
    output = pl.read_csv(result["mixer_input"], separator="\t")

    assert output["N"].to_list() == [1000.0]
    assert result["sample_size_mode"] == "total_n_from_nco"
    assert result["minimum_sample_size"] == 325.0
