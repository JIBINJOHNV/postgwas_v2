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
    FormattingCanonicalColumns,
    FormattingStudyDesign,
    FormattingVcfFields,
)
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.ui.progress import _ACTIVE_PROGRESS
from postgwas.core.variant_identifiers import inspect_bim_identifier_type
from postgwas.modules.formatting.cli import build_parser, main as formatter_main
from postgwas.modules.formatting.contracts import required_formats
from postgwas.modules.formatting.exporters.ldsc import export_ldsc
from postgwas.modules.formatting.exporters.magma import export_magma
from postgwas.modules.formatting.exporters.mixer import export_mixer
from postgwas.modules.formatting.exporters.pred_ld import export_pred_ld
from postgwas.modules.formatting.ldsc_reference import (
    select_ldsc_reference_variants,
)
from postgwas.modules.formatting.service import (
    _configured_schema_reports,
    _ldsc_sample_prevalence_screen_fields,
    _print_contract_summary,
    run_formatter_direct,
)
from postgwas.modules.formatting.table import (
    FormattingError,
    infer_study_design,
    load_harmonised_vcf,
    select_variant_identifiers,
    validate_harmonised_vcf_header,
)
from postgwas.modules.magma.analysis import MagmaReferencePreflight
from postgwas.modules.magma.service import MagmaPipelineResources


FIXTURE = Path(__file__).parent / "fixtures" / "formatting" / "harmonised.vcf"
TARGETS = ["magma", "gcta_gene", "susie", "finemap", "pred_ld", "ldsc", "mixer"]


def _magma_pipeline_context(
    configuration,
    reference_prefix,
    *,
    include_gene_sets=True,
    gene_set_plans=None,
    primary_gene_location_metrics=None,
    analysis_scope=None,
):
    scope = analysis_scope or {
        "mhc_policy": "include",
        "exclude_mhc_snps": False,
        "mhc_region": None,
        "mhc_source": None,
        "exclude_chromosomes": [],
    }
    resources = MagmaPipelineResources(
        configuration=configuration,
        reference=MagmaReferencePreflight(
            ld_reference_prefix=str(reference_prefix),
            executable=sys.executable,
            version="1.10",
            gene_set_plans=gene_set_plans or {},
            analysis_scope=scope,
        ),
        include_gene_sets=include_gene_sets,
        primary_gene_location_metrics=primary_gene_location_metrics,
        file_identities=(),
        log_events=(),
    )
    return RunContext(validations={
        "magma": PipelinePreflightEvidence(
            module="magma",
            input_vcf={},
            resources=resources,
        ),
    })


def test_pipeline_bcftools_resolution_uses_canonical_configuration(
    tmp_path, monkeypatch,
):
    from postgwas.pipeline import runners

    run_config = tmp_path / "run.yaml"
    run_config.write_text(
        "resources:\n"
        "  executables:\n"
        "    bcftools: configured-bcftools\n",
        encoding="utf-8",
    )
    requested = []

    def resolve(command):
        requested.append(command)
        return "/resolved/bin/bcftools"

    monkeypatch.setattr(runners.shutil, "which", resolve)
    args = Namespace(
        run_config=str(run_config),
        bcftools="ignored-former-cli-value",
    )

    runners.check_and_resolve_binaries(args, ["bcftools"])

    assert requested == ["configured-bcftools"]
    assert args.bcftools == "/resolved/bin/bcftools"


def _formatter_header() -> str:
    return FIXTURE.read_text(encoding="utf-8").split("#CHROM", 1)[0]


def _mock_formatter_projection(monkeypatch, values):
    from postgwas.modules.formatting import table

    selected_datasets = []

    monkeypatch.setattr(
        table,
        "read_vcf_header",
        lambda *args, **kwargs: _formatter_header(),
    )

    def write_projection(_vcf, destination, _dataset, columns, _bcftools, **kwargs):
        selected_datasets.append(_dataset)
        names = list(columns)
        Path(destination).write_text(
            "\t".join(names) + "\n"
            + "\t".join(str(values[name]) for name in names) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(table, "extract_vcf_table", write_projection)
    return selected_datasets


def _canonical_projection_values(**overrides):
    values = {
        "CHROM": "1", "POS": "100", "ID": "rs1", "REF": "A", "ALT": "G",
        "BETA": "0.2", "SE": "0.1", "Z": "2", "LP": "8", "EAF": "0.8",
        "INFO": "0.9", "N": "1000", "NEFF": "800", "N_CASE": "400",
        "N_CONTROL": "600",
    }
    values.update(overrides)
    return values


def test_vcf_projection_uses_configured_canonical_names():
    projection = FormattingVcfFields.model_validate({"chromosome_value": "%CHROM"})
    assert projection.root == {"chromosome_value": "%CHROM"}


def test_magma_export_applies_reference_and_scope_before_single_final_write(
    tmp_path,
):
    configuration = load_configuration()
    formatting = configuration.modules.formatting
    magma = configuration.modules.magma
    policy = magma.snp_harmonisation.model_copy(update={
        "resolve_variants_to_reference": True,
    })
    magma = magma.model_copy(update={"snp_harmonisation": policy})
    prefix = tmp_path / "reference"
    (tmp_path / "reference.bim").write_text(
        "1\t1_100_A_G\t0\t100\tA\tG\n"
        "Y\tY_200_C_T\t0\t200\tC\tT\n"
        "6\t6_300_G_A\t0\t300\tG\tA\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "SNP": ["1_100_A_G", "1_150_A_C", "Y_200_C_T", "6_300_G_A"],
        "CHROM": ["1", "1", "Y", "6"],
        "POS": [100, 150, 200, 300],
        "REF": ["A", "A", "C", "G"],
        "ALT": ["G", "C", "T", "A"],
        "LP": [2.0, 3.0, 4.0, 5.0],
        "N": [1000.0, 1000.0, 1000.0, 1000.0],
    })
    result = export_magma(
        frame,
        tmp_path,
        "STUDY",
        formatting,
        overwrite=True,
        magma_config=magma,
        ld_reference_prefix=prefix,
        analysis_scope={
            "mhc_policy": "exclude_both",
            "exclude_mhc_snps": True,
            "mhc_region": {"chromosome": "6", "start": 250, "end": 350},
            "mhc_source": "test",
            "exclude_chromosomes": ["Y", "MT"],
        },
        identifier_qc={
            "identifier_duplicate_groups": 0,
            "identifier_duplicate_rows": 0,
            "identifier_duplicate_rows_excluded": 0,
        },
    )

    locations = pl.read_csv(
        result["snp_loc_file"],
        separator="\t",
        has_header=False,
        new_columns=["SNP", "CHR", "BP"],
    )
    associations = pl.read_csv(result["pval_file"], separator="\t")
    assert locations.to_dicts() == [
        {"SNP": "1_100_A_G", "CHR": 1, "BP": 100},
    ]
    assert associations["SNP"].to_list() == ["1_100_A_G"]
    assert associations["P"].item() == pytest.approx(0.01)
    qc = result["variant_preparation"]["qc"]
    assert qc["reference_unique_id_matches"] == 3
    assert qc["not_in_reference_rows"] == 1
    assert qc["excluded_chromosome_rows"] == 1
    assert qc["excluded_mhc_rows"] == 1
    assert qc["retained_rows"] == 1


def test_magma_export_rejects_duplicate_matching_bim_identifiers(tmp_path):
    configuration = load_configuration()
    formatting = configuration.modules.formatting
    magma = configuration.modules.magma
    prefix = tmp_path / "reference"
    (tmp_path / "reference.bim").write_text(
        "1\t1_100_A_G\t0\t100\tA\tG\n"
        "1\t1_100_A_G\t0\t100\tA\tG\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "SNP": ["1_100_A_G"],
        "CHROM": ["1"],
        "POS": [100],
        "REF": ["A"],
        "ALT": ["G"],
        "LP": [2.0],
        "N": [1000.0],
    })

    with pytest.raises(
        FormattingError,
        match="BIM assigns duplicate records to 1 GWAS variant identifiers",
    ):
        export_magma(
            frame,
            tmp_path,
            "STUDY",
            formatting,
            overwrite=True,
            magma_config=magma,
            ld_reference_prefix=prefix,
            analysis_scope={
                "mhc_policy": "include",
                "exclude_mhc_snps": False,
                "mhc_region": None,
                "mhc_source": None,
                "exclude_chromosomes": [],
            },
        )


def test_pipeline_formatter_checkpoints_magma_reference_preparation(
    tmp_path, monkeypatch,
):
    class CapturingProgress:
        def __init__(self):
            self.started = []

        def start(self, stage):
            self.started.append(stage)

    prefix = tmp_path / "reference"
    bim = tmp_path / "reference.bim"
    bim.write_text("1\t1_100_A_G\t0\t100\tA\tG\n", encoding="utf-8")
    output = tmp_path / "output"
    configuration = load_configuration(cli_overrides={
        "run.dataset_id": "STUDY",
        "run.output_directory": str(output),
        "run.overwrite": True,
        "modules.formatting.formats": ["magma"],
        "modules.formatting.variant_identifiers.target_types": {
            "magma": "unique",
        },
        "modules.formatting.variant_identifiers.target_duplicate_policies": {
            "magma": "lowest_p",
        },
        "modules.magma.input.ld_reference_prefix": str(prefix),
    })
    frame = pl.DataFrame({
        "CHROM": ["1"],
        "POS": [100],
        "ID": ["rs1"],
        "REF": ["A"],
        "ALT": ["G"],
        "LP": [2.0],
        "N": [1000.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )
    args = Namespace(
        vcf=str(FIXTURE),
        bcftools=sys.executable,
        _pipeline_progress_plan={
            "kind": "pops",
            "magma_stage_numbers": {"variant_inputs": 4},
        },
        _pipeline_stage_progress=CapturingProgress(),
    )
    analysis_scope = {
            "mhc_policy": "include",
            "exclude_mhc_snps": False,
            "mhc_region": None,
            "mhc_source": None,
            "exclude_chromosomes": [],
    }

    result = run_formatter_direct(
        args,
        ctx=_magma_pipeline_context(
            configuration,
            prefix,
            include_gene_sets=False,
            analysis_scope=analysis_scope,
        ),
        configuration=configuration,
        emit_terminal_summary=False,
    )["magma"]

    manifest = yaml.safe_load(
        (output / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8",
        )
    )
    assert manifest["input_resources"]["magma_bim"]["path"] == str(
        bim.resolve()
    )
    assert result["variant_preparation"]["qc"]["retained_rows"] == 1
    assert Path(result["variant_preparation"]["excluded_variants"]).is_file()
    assert args._pipeline_stage_progress.started == [4]

    resumed_configuration = load_configuration(cli_overrides={
        "run.dataset_id": "STUDY",
        "run.output_directory": str(output),
        "modules.formatting.formats": ["magma"],
        "modules.formatting.variant_identifiers.target_types": {
            "magma": "unique",
        },
        "modules.formatting.variant_identifiers.target_duplicate_policies": {
            "magma": "lowest_p",
        },
        "modules.magma.input.ld_reference_prefix": str(prefix),
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail(
            "validated MAGMA formatter outputs should resume without VCF extraction"
        ),
    )
    resumed_args = Namespace(
        vcf=str(FIXTURE),
        bcftools=sys.executable,
        _pipeline_progress_plan=args._pipeline_progress_plan,
        _pipeline_stage_progress=CapturingProgress(),
    )
    resumed = run_formatter_direct(
        resumed_args,
        ctx=_magma_pipeline_context(
            resumed_configuration,
            prefix,
            include_gene_sets=False,
            analysis_scope=analysis_scope,
        ),
        configuration=resumed_configuration,
        emit_terminal_summary=False,
    )["magma"]
    assert resumed["variant_preparation"]["qc"] == result[
        "variant_preparation"
    ]["qc"]
    assert resumed_args._pipeline_stage_progress.started == [4]


def test_formatter_input_contract_matches_harmonisation_output_headers():
    configuration = load_configuration()
    formatting = configuration.modules.formatting.input_contract
    harmonisation = configuration.modules.harmonisation.vcf_processing

    assert (
        "##" + formatting.genome_build_metadata
        == harmonisation.genome_build_header
    )
    produced = harmonisation.provenance.headers
    assert formatting.provenance_headers.model_dump() == {
        "version": produced["postgwas_version"],
        "dataset_id": produced["dataset_id"],
        "status": produced["vcf_status"],
    }


def test_formatter_header_requires_postgwas_provenance():
    config = load_configuration().modules.formatting
    header = _formatter_header().replace('##postgwas_version="test"\n', "")

    with pytest.raises(
        FormattingError,
        match=r"missing PostGWAS provenance header\(s\): postgwas_version",
    ):
        validate_harmonised_vcf_header(header, config)


def test_formatter_rejects_invalid_header_before_extracting_records(
    tmp_path, monkeypatch,
):
    from postgwas.modules.formatting import table

    config = load_configuration().modules.formatting
    invalid_header = _formatter_header().replace(
        '##postgwas_vcf_status="created and structurally validated"\n',
        "",
    )
    monkeypatch.setattr(
        table,
        "read_vcf_header",
        lambda *args, **kwargs: invalid_header,
    )
    monkeypatch.setattr(
        table,
        "extract_vcf_table",
        lambda *args, **kwargs: pytest.fail(
            "formatter extracted records before validating the VCF header"
        ),
    )

    with pytest.raises(
        FormattingError,
        match=r"missing PostGWAS provenance header\(s\): postgwas_vcf_status",
    ):
        load_harmonised_vcf(
            tmp_path / "input.vcf",
            tmp_path / "projection.tsv",
            "bcftools",
            config,
        )


def test_formatter_records_declared_dataset_without_requiring_output_id_match():
    config = load_configuration().modules.formatting

    evidence = validate_harmonised_vcf_header(_formatter_header(), config)

    assert evidence["postgwas_dataset_id"] == "STUDY"


def test_formatter_load_uses_vcf_declared_dataset_for_sample_selection(
    tmp_path, monkeypatch,
):
    config = load_configuration().modules.formatting
    selected_datasets = _mock_formatter_projection(
        monkeypatch, _canonical_projection_values(),
    )

    frame = load_harmonised_vcf(
        tmp_path / "input.vcf",
        tmp_path / "projection.tsv",
        "bcftools",
        config,
    )

    assert frame.height == 1
    assert selected_datasets == ["STUDY"]


def test_formatter_load_can_return_validated_header_evidence(
    tmp_path, monkeypatch,
):
    config = load_configuration().modules.formatting
    _mock_formatter_projection(monkeypatch, _canonical_projection_values())

    frame, evidence = load_harmonised_vcf(
        tmp_path / "input.vcf",
        tmp_path / "projection.tsv",
        "bcftools",
        config,
        return_header_evidence=True,
    )

    assert frame.height == 1
    assert evidence["genome_build"] == "GRCh37"
    assert evidence["postgwas_dataset_id"] == "STUDY"
    assert evidence["postgwas_version"] == "test"
    assert evidence["postgwas_status"] == "created and structurally validated"


def test_formatter_reuses_pipeline_validated_header_evidence(
    tmp_path, monkeypatch,
):
    from postgwas.modules.formatting import table

    config = load_configuration().modules.formatting
    selected_datasets = []
    evidence = validate_harmonised_vcf_header(_formatter_header(), config)

    monkeypatch.setattr(
        table,
        "read_vcf_header",
        lambda *args, **kwargs: pytest.fail(
            "pipeline-validated VCF header was read a second time"
        ),
    )

    def write_projection(_vcf, destination, dataset, columns, _bcftools, **kwargs):
        selected_datasets.append(dataset)
        values = _canonical_projection_values()
        Path(destination).write_text(
            "\t".join(columns) + "\n"
            + "\t".join(str(values[name]) for name in columns) + "\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(table, "extract_vcf_table", write_projection)
    frame, returned = load_harmonised_vcf(
        tmp_path / "input.vcf",
        tmp_path / "projection.tsv",
        "bcftools",
        config,
        return_header_evidence=True,
        validated_header_evidence=evidence,
    )

    assert frame.height == 1
    assert returned == evidence
    assert selected_datasets == ["STUDY"]


def test_formatter_parses_numeric_text_without_rewriting_canonical_values(
    tmp_path, monkeypatch,
):
    config = load_configuration().modules.formatting
    _mock_formatter_projection(monkeypatch, _canonical_projection_values())

    frame = load_harmonised_vcf(
        tmp_path / "input.vcf",
        tmp_path / "projection.tsv",
        "bcftools",
        config,
    )

    assert frame["CHROM"].to_list() == ["1"]
    assert frame["REF"].to_list() == ["A"]
    assert frame["ALT"].to_list() == ["G"]
    assert frame.schema["POS"] == pl.Int64
    assert frame["POS"].to_list() == [100]
    assert all(
        frame.schema[column] == pl.Float64
        for column in config.numeric_columns
    )


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({"CHROM": "chr1"}, "1 noncanonical chromosome labels"),
        ({"REF": "a"}, "1 noncanonical or invalid REF/ALT allele records"),
        ({"POS": "not-a-position"}, "1 invalid positions"),
    ],
)
def test_formatter_rejects_noncanonical_rows_instead_of_repairing_them(
    tmp_path, monkeypatch, overrides, expected,
):
    config = load_configuration().modules.formatting
    _mock_formatter_projection(
        monkeypatch,
        _canonical_projection_values(**overrides),
    )

    with pytest.raises(FormattingError, match=expected):
        load_harmonised_vcf(
            tmp_path / "input.vcf",
            tmp_path / "projection.tsv",
            "bcftools",
            config,
        )


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
    assert "Trait-type inference" in screen_text
    assert "not needed to prepare the selected input files" in screen_text
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

    assert _ACTIVE_PROGRESS.get() == ()


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
    assert result["ldsc"]["sample_prevalence_aggregation"] == "median"
    assert result["ldsc"]["sample_prevalence_variants"] == 3
    assert result["ldsc"]["sample_prevalence_minimum"] == pytest.approx(
        350 / (350 + 550)
    )
    assert result["ldsc"]["sample_prevalence_maximum"] == pytest.approx(0.4)
    assert result["ldsc"]["sample_prevalence_case_count_minimum"] == 350
    assert result["ldsc"]["sample_prevalence_case_count_maximum"] == 400
    assert result["ldsc"]["sample_prevalence_control_count_minimum"] == 550
    assert result["ldsc"]["sample_prevalence_control_count_maximum"] == 600

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
    assert resolved["modules"]["formatting"]["ldsc_sample_prevalence"] == {
        "aggregation": "median",
    }
    manifest = yaml.safe_load(
        (tmp_path / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    manifest_ldsc = manifest["results"]["ldsc"]
    assert manifest_ldsc["sample_prev"] == pytest.approx(0.4)
    assert manifest_ldsc["sample_prevalence_aggregation"] == "median"
    assert manifest_ldsc["sample_prevalence_variants"] == 3
    assert manifest_ldsc["sample_prevalence_case_count_minimum"] == 350
    assert manifest_ldsc["sample_prevalence_case_count_maximum"] == 400
    assert manifest_ldsc["sample_prevalence_control_count_minimum"] == 550
    assert manifest_ldsc["sample_prevalence_control_count_maximum"] == 600
    log_text = (tmp_path / "logs" / "STUDY_formatter.log").read_text(encoding="utf-8")
    assert "formatter_run status=COMPLETED" in log_text
    assert "formatter_input_vcf_contract status=PASSED" in log_text
    assert "genome_build=GRCh37" in log_text
    assert "postgwas_dataset_id=STUDY" in log_text
    assert "postgwas_version=test" in log_text
    assert "formatter_input_row_contract status=PASSED" in log_text
    assert "values_rewritten=0" in log_text
    assert "formatter_saved_schema target=ldsc" in log_text
    assert "target=ldsc variant_id_type=rsid" in log_text
    assert "p_value_saved='LP → P (raw P (10^(-LP)))'" in log_text
    assert "allele_frequency_saved='EAF → FRQ" in log_text
    assert "(effect-allele frequency (EAF))'" in log_text
    assert "sample_size_saved='case N: FORMAT/NC → N_CAS (copied);" in log_text
    assert "control N: FORMAT/NCO → N_CON" in log_text
    assert "total/effective N are not written'" in log_text
    assert "DECIDE   study_design trait_type=binary" in log_text
    assert "DECIDE   ldsc_sample_prevalence aggregation=median" in log_text
    assert "formula=N_CASE/(N_CASE+N_CONTROL)" in log_text
    assert "variants_used=3" in log_text
    assert "case_count_minimum=350" in log_text
    assert "case_count_maximum=400" in log_text
    assert "control_count_minimum=550" in log_text
    assert "control_count_maximum=600" in log_text
    assert "metadata_header_used=false" in log_text
    screen_text = capsys.readouterr().out
    normalized_screen = " ".join(screen_text.split())
    assert "Formatter plan · prepare files required by downstream analyses" in screen_text
    assert "GCTA COJO / fastBAT / mBAT-combo" in screen_text
    assert "Variant identifiers written" in normalized_screen
    assert "rsID" in normalized_screen
    assert "Columns and transformations" in normalized_screen
    assert "LP → P" in normalized_screen
    assert "EAF → maf" in normalized_screen
    assert "N_CASE" in normalized_screen
    assert "case N: FORMAT/NC" in normalized_screen
    assert "N_CAS" in normalized_screen
    assert "N_CONTROL" in normalized_screen
    assert "control N:" in normalized_screen
    assert "FORMAT/NCO" in normalized_screen
    assert "N_CON" in normalized_screen
    assert "copied" in normalized_screen
    assert "minor-allele" in normalized_screen
    assert "effect-allele" in normalized_screen
    assert "Source GWAS-VCF" in screen_text
    assert "Validate and read the harmonised GWAS-VCF" in screen_text
    assert "Total variants in input GWAS-VCF" in screen_text
    assert "Total variants in input GWAS-VCF : 4" in normalized_screen
    assert "Genome build" in screen_text
    assert "Genome build : GRCh37" in normalized_screen
    assert "VCF embedded dataset/sample" in screen_text
    assert "VCF embedded dataset/sample : STUDY" in normalized_screen
    assert "Run ID versus VCF identity" in screen_text
    assert "Run ID versus VCF identity : match" in normalized_screen
    assert "Trait type" in screen_text
    assert "Binary case-control (inferred from sample-count fields)" in screen_text
    assert "N_CASE completeness" in screen_text
    assert "4 / 4 variants (100%)" in screen_text
    assert "LDSC files created" in screen_text
    assert "Variants written to LDSC files" in screen_text
    assert "Required-value validation" in screen_text
    assert "Interpretation" in screen_text
    assert "downstream analysis can retain fewer variants" in normalized_screen
    assert "Prepared scientific files" in screen_text
    assert Path(result["magma"]["snp_loc_file"]).name in screen_text
    assert Path(result["magma"]["pval_file"]).name in screen_text
    assert "GWAS-VCF case fraction" in screen_text
    assert "Calculated case fraction" in screen_text
    assert "0.4 (40.00%)" in screen_text
    assert "N_CASE / (N_CASE + N_CONTROL)" in screen_text
    assert "median across retained LDSC-input variants" in screen_text
    assert "Variants used" in screen_text
    assert "0.3888888889 to 0.4" in screen_text
    assert "N_CASE range" in screen_text
    assert "350 to 400" in screen_text
    assert "N_CONTROL range" in screen_text
    assert "550 to 600" in screen_text
    assert "Detailed HTML report" in screen_text
    screen_lines = screen_text.splitlines()

    def summary_line(label):
        return next(line for line in screen_lines if label in line)

    assert summary_line("input files prepared").startswith("  🔬")
    assert summary_line("Source GWAS-VCF").startswith("      🔬")
    assert summary_line("Dataset").startswith("          🔹")
    assert summary_line("LDSC files created").startswith("      🔬")
    assert summary_line("Variant accounting").startswith(
        "          🔬"
    )
    assert summary_line("Variants read from GWAS-VCF").startswith("              🧮")
    assert summary_line("GWAS-VCF case fraction").startswith("          🔬")
    assert summary_line("Calculated case fraction").startswith("              🔬")
    assert summary_line("Reports and logs").startswith("      🔬")
    assert summary_line("Detailed HTML report").startswith("          ✅")
    aligned_fields = [
        summary_line(label)
        for label in (
            "Dataset",
            "Input GWAS-VCF",
            "Variants read from GWAS-VCF",
            "Variants written to LDSC files",
            "Final identifier validation",
            "Detailed HTML report",
            "Full formatter log",
        )
    ]
    assert len({line.index(" : ") for line in aligned_fields}) == 1
    html_report = Path(result["ldsc"]["html_report"])
    assert html_report == tmp_path / "reports" / "STUDY_formatter_report.html"
    assert html_report.is_file() and html_report.stat().st_size > 0
    assert all(
        Path(target_result["html_report"]) == html_report
        for target_result in result.values()
    )
    html = html_report.read_text(encoding="utf-8")
    assert "Formatter report" in html
    assert "Target export overview" in html
    assert "Input and study-design evidence" in html
    assert "Identifier and duplicate handling" in html
    assert "Configured output schema" in html
    assert "Scientific statistic semantics" in html
    assert "Complete artifact inventory" in html
    assert "Execution provenance" in html
    assert "Evidence reuse" in html
    assert "STUDY_ldsc_input.tsv" in html
    assert "href=\"../STUDY_ldsc_input.tsv\"" in html
    assert str(tmp_path) not in html
    assert any(
        record["path"] == str(html_report)
        for record in manifest["outputs"]
    )
    summary_lines = [
        line
        for line in screen_text.splitlines()
        if any(
            label in line
            for label in (
                "Dataset",
                "Input GWAS-VCF",
                "Trait type",
                "Detailed HTML report",
                "Full formatter log",
            )
        )
        and ":" in line
    ]
    assert len(summary_lines) == 5
    assert len({line.index(":") for line in summary_lines}) == 1


def test_ldsc_sample_prevalence_is_printed_after_validated_resume(
    tmp_path, monkeypatch, capsys,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "1"],
        "POS": [100, 200, 300],
        "ID": ["rs1", "rs2", "rs3"],
        "REF": ["A", "C", "G"],
        "ALT": ["G", "T", "A"],
        "Z": [2.0, -2.0, 3.0],
        "LP": [8.0, 6.0, 10.0],
        "N_CASE": [100.0, 400.0, 900.0],
        "N_CONTROL": [900.0, 600.0, 100.0],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )
    common = {
        "vcf": str(FIXTURE),
        "dataset_id": "STUDY",
        "output_directory": str(tmp_path),
        "format": ["ldsc"],
        "bcftools": sys.executable,
        "resume": True,
    }

    run_formatter_direct(Namespace(**common, overwrite=True))
    capsys.readouterr()
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: pytest.fail("resume re-read the GWAS-VCF"),
    )

    resumed = run_formatter_direct(Namespace(**common))
    screen_text = " ".join(capsys.readouterr().out.split())

    assert resumed["ldsc"]["resumed"] is True
    assert "Formatter outputs validated" in screen_text
    assert "Variant identifiers written" in screen_text
    assert "rsID" in screen_text
    assert "GWAS-VCF case fraction" in screen_text
    assert "Calculated case fraction" in screen_text
    assert "0.4 (40.00%)" in screen_text
    assert "N_CASE / (N_CASE + N_CONTROL)" in screen_text
    assert "median across retained LDSC-input variants" in screen_text
    assert "Variants used" in screen_text
    assert "0.1 to 0.9" in screen_text
    assert "N_CASE range" in screen_text
    assert "100 to 900" in screen_text
    assert "N_CONTROL range" in screen_text
    assert "100 to 900" in screen_text


def test_legacy_ldsc_resume_metadata_requests_rerun_for_count_ranges():
    fields = _ldsc_sample_prevalence_screen_fields(
        {
            "ldsc": {
                "trait_type": "binary",
                "sample_prev": 0.25,
                "sample_prevalence_aggregation": "median",
                "sample_prevalence_variants": 100,
                "sample_prevalence_minimum": 0.2,
                "sample_prevalence_maximum": 0.3,
            },
        },
        "N_CASE",
        "N_CONTROL",
        48,
    )
    screen_text = " ".join(fields)

    assert "Calculated case fraction" in screen_text
    assert "0.25 (25.00%)" in screen_text
    assert "N_CASE/N_CONTROL values" in screen_text
    assert "rerun with --overwrite" in screen_text


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
    current_report = (
        copied / "reports" / "STUDY_formatter_report.html"
    ).resolve()
    result = resumed["gcta_gene"]
    assert Path(result["summary_statistics_input_file"]) == current_artifact
    assert Path(result["outputs"]["summary_statistics"]["path"]) == current_artifact
    assert Path(result["log_file"]) == current_log
    assert Path(result["html_report"]) == current_report
    assert all(
        Path(target_result["html_report"]) == current_report
        for target_result in resumed.values()
    )
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
    copied_html = current_report.read_text(encoding="utf-8")
    assert str(original) not in copied_html
    assert "href=\"../STUDY_gcta.ma\"" in copied_html

    copied_manifest = yaml.safe_load(
        (copied / "run_metadata" / "formatter_completion.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert copied_manifest["results"]["gcta_gene"][
        "summary_statistics_input_file"
    ] == str(current_artifact)
    assert copied_manifest["results"]["gcta_gene"]["log_file"] == str(current_log)
    assert copied_manifest["results"]["gcta_gene"]["html_report"] == str(
        current_report
    )
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
    assert {record["path"] for record in manifest["outputs"]} == {
        str(current_artifact),
        resumed["custom"]["html_report"],
    }


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
def test_formatter_resume_restarts_when_all_current_artifacts_are_missing(tmp_path):
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

    restarted = run_formatter_direct(Namespace(
        **common, output_directory=str(copied),
    ))

    assert original_artifact.is_file()
    assert copied_artifact.is_file()
    assert restarted["gcta_gene"].get("resumed") is not True
    assert Path(
        restarted["gcta_gene"]["summary_statistics_input_file"]
    ) == copied_artifact
    log = copied / "logs" / "STUDY_formatter.log"
    assert "reason=incomplete_outputs" in log.read_text(
        encoding="utf-8"
    )


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
def test_formatter_resume_restarts_changed_input_vcf(tmp_path):
    vcf = tmp_path / "input.vcf"
    shutil.copyfile(FIXTURE, vcf)
    output = tmp_path / "output"
    common = dict(
        vcf=str(vcf), dataset_id="STUDY", output_directory=str(output),
        format=["gcta_gene"], bcftools=shutil.which("bcftools"), resume=True,
    )
    run_formatter_direct(Namespace(**common, overwrite=True))
    original = vcf.read_text(encoding="utf-8")
    changed = original.replace(
        "0.2:0.1:2:8:0.8:0.9:1000:800:400:600",
        "0.25:0.1:2:8:0.8:0.9:1000:800:400:600",
        1,
    )
    assert changed != original
    vcf.write_text(changed, encoding="utf-8")

    restarted = run_formatter_direct(Namespace(**common))

    assert restarted["gcta_gene"].get("resumed") is not True
    log = output / "logs" / "STUDY_formatter.log"
    log_text = log.read_text(encoding="utf-8")
    assert "formatter input VCF changed since this checkpoint" in log_text
    assert "reason=changed_inputs" in log_text


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_formatter_resume_restarts_changed_parameters(tmp_path):
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["gcta_gene"],
        bcftools=shutil.which("bcftools"),
        resume=True,
    )
    run_formatter_direct(Namespace(**common, overwrite=True))

    restarted = run_formatter_direct(Namespace(
        **common,
        variant_id_types={"gcta_gene": "unique"},
    ))

    assert restarted["gcta_gene"].get("resumed") is not True
    log_text = (
        tmp_path / "logs" / "STUDY_formatter.log"
    ).read_text(encoding="utf-8")
    assert "Resolved formatter parameters changed" in log_text
    assert "reason=changed_parameters" in log_text


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
@pytest.mark.parametrize("target", ["gcta_gene", "magma"])
def test_formatter_refuses_outputs_without_checksum_manifest(
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
    with pytest.raises(
        FormattingError,
        match="no checksum-validated completion manifest",
    ):
        run_formatter_direct(Namespace(**common))

    assert not manifest.exists()


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

    with pytest.raises(
        FormattingError,
        match="no checksum-validated completion manifest",
    ):
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
        (
            "reports/STUDY_formatter_report.html",
            "formatter.html_report and custom",
        ),
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

    expected = 1 + sum(
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
    assert destinations["formatter.html_report"] == str(
        tmp_path / "reports" / "STUDY_formatter_report.html"
    )


def test_formatter_html_report_path_is_schema_validated():
    runtime = load_configuration().modules.formatting.runtime

    assert runtime.html_report_file == "reports/{dataset_id}_formatter_report.html"
    invalid = runtime.model_dump(mode="python")
    invalid["html_report_file"] = "reports/{dataset_id}_formatter_report.txt"
    with pytest.raises(ValueError, match="html_report_file must end with .html"):
        runtime.__class__.model_validate(invalid)


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
    assert "Outside configured chromosomes" in screen_text
    assert "1 (50.00% of GWAS-VCF variants)" in screen_text
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
    normalized_help = " ".join(help_text.split())
    output_action = next(
        action
        for action in build_parser()._actions
        if action.dest == "output_directory"
    )
    assert "The VCF is read once" in help_text
    assert "--vcf PATH --output-directory PATH" in normalized_help
    assert output_action.required is True
    assert "pred_ld" in help_text
    assert "mixer" in help_text
    assert "modules.formatting.formats" in help_text
    assert "--variant-id-type {rsid,unique}" in help_text
    assert "--merge-alleles PATH" in help_text
    assert "direct formatter uses it optionally" in normalized_help
    assert "selected duplicate-ID policy handles duplicated rsIDs" in normalized_help
    assert "Per-target YAML settings" in help_text
    assert "--run-config PATH" in help_text
    assert "--custom-output FILE" in help_text
    assert "--id NAME" in help_text
    assert "--n-case NAME" in help_text
    assert "without editing YAML" in help_text


def test_every_formatter_target_reports_exact_columns_p_and_frequency_semantics():
    module = load_configuration().modules.formatting
    reports = _configured_schema_reports(
        TARGETS,
        module,
        {"ldsc": "binary", "mixer": "binary"},
    )

    assert list(reports) == TARGETS
    assert all(
        report["variant_id_type"] == "rsid"
        for report in reports.values()
    )
    assert all(
        report["variant_id_type_label"] == "rsID"
        for report in reports.values()
    )
    assert reports["magma"]["column_mappings"].endswith(
        "p values: SNP → SNP, LP → P, N → N_COL"
    )
    assert reports["magma"]["p_value"] == (
        "p values: LP → P (raw P (10^(-LP)))"
    )
    assert reports["gcta_gene"]["frequency"] == (
        "EAF → freq (effect-allele frequency (EAF))"
    )
    assert reports["gcta_gene"]["name"] == (
        "GCTA COJO / fastBAT / mBAT-combo"
    )
    assert reports["magma"]["sample_size"] == (
        "p values · total N: FORMAT/SS → N_COL (copied); "
        "used as MAGMA per-variant total N"
    )
    assert reports["gcta_gene"]["sample_size"] == (
        "total N: FORMAT/SS → N (copied); used as GCTA per-variant total N"
    )
    assert reports["susie"]["p_value"] == "LP → LP (-log10(P))"
    assert reports["susie"]["sample_size"] == (
        "harmonised N (binary effective N; quantitative total N): "
        "FORMAT/NEF → NEF (copied); fine-mapping uses the locus median as "
        "SuSiE n"
    )
    assert reports["finemap"]["frequency"] == (
        "EAF → maf (minor-allele frequency (MAF = min(EAF, 1-EAF)))"
    )
    assert reports["finemap"]["sample_size"] == (
        "harmonised N (binary effective N; quantitative total N): "
        "FORMAT/NEF → NEF (copied); fine-mapping uses the rounded locus "
        "median as FINEMAP n_samples"
    )
    assert reports["pred_ld"]["p_value"] == "LP → LP (-log10(P))"
    assert reports["pred_ld"]["frequency"] == (
        "EAF → AF (effect-allele frequency (EAF))"
    )
    assert reports["pred_ld"]["sample_size"] == (
        "binary control N / quantitative total N: FORMAT/NCO → NC (copied); "
        "total N: FORMAT/SS → SS (copied); not consumed by PRED-LD; retained "
        "for re-harmonisation"
    )
    assert reports["ldsc"]["p_value"] == "LP → P (raw P (10^(-LP)))"
    assert reports["ldsc"]["frequency"] == (
        "EAF → FRQ (effect-allele frequency (EAF))"
    )
    assert reports["ldsc"]["sample_size"] == (
        "case N: FORMAT/NC → N_CAS (copied); control N: FORMAT/NCO → N_CON "
        "(copied); case/control counts are written separately; total/effective "
        "N are not written"
    )
    assert "N_CASE → N_CAS" in reports["ldsc"]["column_mappings"]
    assert reports["mixer"]["p_value"] == "not saved"
    assert reports["mixer"]["frequency"] == "not saved"
    assert reports["mixer"]["sample_size"] == (
        "effective N = 4 / (1 / cases + 1 / controls): FORMAT/NEF → N "
        "(copied); MiXeR uses saved N and excludes variants below the "
        "configured fraction of median N"
    )


def test_quantitative_ldsc_and_custom_schema_reports_use_resolved_names():
    module = load_configuration().modules.formatting
    variant_identifiers = module.variant_identifiers.model_copy(update={
        "target_types": {"ldsc": "unique"},
    })
    custom = module.custom_output.model_copy(update={
        "output_file": "custom.tsv",
        "columns": {
            "id": "MARKER",
            "p": "PVALUE",
            "lp": "LOGP",
            "eaf": "EFFECT_FREQ",
            "maf": "MINOR_FREQ",
            "n": "TOTAL_N",
            "neff": "EFFECTIVE_N",
            "n_case": "CASES",
            "n_control": "CONTROLS",
        },
    })
    module = module.model_copy(update={
        "custom_output": custom,
        "variant_identifiers": variant_identifiers,
    })

    reports = _configured_schema_reports(
        ["ldsc", "custom"],
        module,
        {"ldsc": "quantitative"},
    )

    assert "N_CONTROL → N" in reports["ldsc"]["column_mappings"]
    assert "N_CASE" not in reports["ldsc"]["column_mappings"]
    assert reports["ldsc"]["variant_id_type"] == "unique"
    assert reports["ldsc"]["variant_id_type_label"] == (
        "coordinate-and-allele ID (CHROM_POS_REF_ALT)"
    )
    assert reports["custom"]["variant_id_type"] == "rsid"
    assert reports["custom"]["variant_id_type_label"] == "rsID"
    assert reports["custom"]["column_mappings"] == (
        "SNP → MARKER, LP → PVALUE, LP → LOGP, EAF → EFFECT_FREQ, "
        "EAF → MINOR_FREQ, N → TOTAL_N, NEFF → EFFECTIVE_N, "
        "N_CASE → CASES, N_CONTROL → CONTROLS"
    )
    assert reports["ldsc"]["sample_size"] == (
        "total N: FORMAT/NCO → N (copied); total N is written; "
        "case/control/effective N are not written"
    )
    assert reports["custom"]["p_value"] == (
        "LP → PVALUE (raw P (10^(-LP))); LP → LOGP (-log10(P))"
    )
    assert reports["custom"]["frequency"] == (
        "EAF → EFFECT_FREQ (effect-allele frequency (EAF)); "
        "EAF → MINOR_FREQ (minor-allele frequency "
        "(MAF = min(EAF, 1-EAF)))"
    )
    assert reports["custom"]["sample_size"] == (
        "total N: FORMAT/SS → TOTAL_N (copied); harmonised N (binary effective "
        "N; quantitative total N): FORMAT/NEF → EFFECTIVE_N (copied); case N: "
        "FORMAT/NC → CASES (copied); binary control N / quantitative total N: "
        "FORMAT/NCO → CONTROLS (copied)"
    )


def test_sample_size_reporting_contract_is_complete_and_trait_aware():
    reporting = load_configuration().modules.formatting.sample_size_reporting

    assert reporting.source_semantics["effective_sample_size"].resolve(
        "binary"
    ) == "effective N = 4 / (1 / cases + 1 / controls)"
    assert reporting.source_semantics["effective_sample_size"].resolve(
        "quantitative"
    ) == "total N (NEF = NCO)"
    assert reporting.target_notes["pred_ld"].default == (
        "not consumed by PRED-LD; retained for re-harmonisation"
    )

    incomplete = reporting.model_dump(mode="python")
    incomplete["source_semantics"].pop("effective_sample_size")
    with pytest.raises(ValueError, match="every sample-size role"):
        reporting.__class__.model_validate(incomplete)


def test_sample_size_roles_require_distinct_canonical_columns():
    canonical = (
        load_configuration().modules.formatting.canonical_columns.model_dump()
    )
    canonical["effective_sample_size"] = canonical["total_sample_size"]

    with pytest.raises(ValueError, match="sample-size roles must reference distinct"):
        FormattingCanonicalColumns.model_validate(canonical)


def test_formatter_contract_summary_is_aligned_and_reports_sample_size(capsys):
    module = load_configuration().modules.formatting
    reports = _configured_schema_reports(
        ["ldsc", "mixer"],
        module,
        {"ldsc": "binary", "mixer": "binary"},
    )

    _print_contract_summary(reports, 42)

    output = capsys.readouterr().out
    screen = " ".join(output.split())
    for expected in (
        "case N:",
        "FORMAT/NC",
        "N_CAS",
        "control N:",
        "FORMAT/NCO",
        "N_CON",
        "effective N = 4 / (1 /",
        "cases + 1 / controls):",
        "FORMAT/NEF",
    ):
        assert expected in screen
    field_lines = [
        line
        for line in output.splitlines()
        if any(
            label in line
            for label in (
                "Variant identifiers written",
                "Duplicate-ID handling",
                "Columns and transformations",
                "P value written",
                "Allele frequency written",
                "Sample size written",
            )
        )
        and " : " in line
    ]
    assert field_lines
    assert len({line.index(" : ") for line in field_lines}) == 1


def test_sample_size_reporting_does_not_change_scientific_checkpoint_digest(
    tmp_path,
):
    from postgwas.modules.formatting.resume import (
        _configuration_digest,
        formatter_resolved_paths,
    )

    baseline = load_configuration()
    override = tmp_path / "reporting_only.yaml"
    override.write_text(
        "modules:\n"
        "  formatting:\n"
        "    sample_size_reporting:\n"
        "      target_notes:\n"
        "        magma:\n"
        "          default: shown with alternative screen wording\n",
        encoding="utf-8",
    )
    changed = load_configuration(override)

    assert "sample_size_reporting" in formatter_resolved_paths(
        changed, ["magma"]
    )
    assert _configuration_digest(baseline, ["magma"]) == (
        _configuration_digest(changed, ["magma"])
    )


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


def test_formatter_requires_explicit_output_directory(capsys):
    with pytest.raises(SystemExit) as caught:
        formatter_main([
            "--vcf", str(FIXTURE),
            "--dataset-id", "STUDY",
            "--format", "ldsc",
        ])

    assert caught.value.code == 2
    assert "--output-directory" in capsys.readouterr().err


def test_formatter_configurable_cli_actions_do_not_own_defaults():
    parser = build_parser()
    for destination in (
        "format", "run_config", "resume", "overwrite", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed", "variant_id_type",
        "duplicate_id_policy", "merge_alleles", "custom_output_file",
        "custom_columns",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS
    assert "bcftools" not in {action.dest for action in parser._actions}


def test_formatter_merge_alleles_is_optional_and_validated(tmp_path):
    without_reference = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--output-directory", str(tmp_path / "without_reference"),
        "--format", "ldsc",
    ])
    assert not hasattr(without_reference, "merge_alleles")

    reference = tmp_path / "w_hm3.snplist"
    reference.write_text("SNP A1 A2\nrs1 G A\n", encoding="utf-8")
    with_reference = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--output-directory", str(tmp_path / "with_reference"),
        "--format", "ldsc",
        "--merge-alleles", str(reference),
    ])
    assert Path(with_reference.merge_alleles) == reference


def test_heritability_requires_ldsc_resources_in_direct_and_pipeline_modes():
    from postgwas.modules.ldsc.cli import build_parser as build_ldsc_parser
    from postgwas.pipeline.registry import REGISTRY

    parser = build_ldsc_parser()
    for destination in (
        "ldsc_input", "merge_alleles", "ref_ld_chr", "w_ld_chr",
    ):
        action = next(
            item for item in parser._actions if item.dest == destination
        )
        assert action.required is True

    required = {
        option.flag
        for option in REGISTRY.get("heritability").required_options
    }
    assert {"--merge-alleles", "--ref-ld-chr", "--w-ld-chr"} <= required


def test_custom_cli_preserves_requested_column_order(tmp_path):
    args = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--output-directory", str(tmp_path / "output"),
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
    assert {record["path"] for record in manifest["outputs"]} == {
        result["custom_output_file"],
        result["html_report"],
    }


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


def test_identifier_selection_can_exclude_every_duplicated_id_row():
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "1", "2"],
        "POS": [100, 101, 200, 300],
        "ID": ["rs1", "rs1", "rs2", "."],
        "REF": ["A", "A", "C", "G"],
        "ALT": ["G", "C", "T", "A"],
    })
    module = load_configuration().modules.formatting

    selected, qc = select_variant_identifiers(
        frame,
        module,
        "rsid",
        duplicate_policy="exclude_all",
    )

    assert selected["SNP"].to_list() == ["rs2"]
    assert qc == {
        "variant_id_type": "rsid",
        "identifier_rows_in": 4,
        "identifier_rows_out": 1,
        "identifier_rows_excluded": 3,
        "identifier_missing_rows_excluded": 1,
        "identifier_duplicate_groups": 1,
        "identifier_duplicate_rows": 2,
        "identifier_duplicate_policy": "exclude_all",
        "identifier_exact_duplicate_rows_collapsed": 0,
        "identifier_conflicting_duplicate_groups": 1,
        "identifier_conflicting_duplicate_rows": 2,
        "identifier_duplicate_groups_resolved_by_policy": 0,
        "identifier_duplicate_groups_unresolved": 1,
        "identifier_ranked_winners_retained": 0,
        "identifier_duplicate_rows_excluded": 2,
    }


def test_identifier_selection_rejects_when_all_ids_are_duplicated():
    frame = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 101],
        "ID": ["rs1", "rs1"],
        "REF": ["A", "A"],
        "ALT": ["G", "C"],
    })
    module = load_configuration().modules.formatting

    with pytest.raises(
        FormattingError,
        match="No variants remain.*Provide a compatible reference",
    ):
        select_variant_identifiers(
            frame,
            module,
            "rsid",
            duplicate_policy="exclude_all",
        )


def test_identifier_selection_collapses_only_exact_repeated_records():
    frame = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 100],
        "ID": ["rs1", "rs1"],
        "REF": ["A", "A"],
        "ALT": ["G", "G"],
        "LP": [8.0, 8.0],
    })
    module = load_configuration().modules.formatting

    selected, qc = select_variant_identifiers(
        frame, module, "rsid", duplicate_policy="error",
    )

    assert selected.height == 1
    assert qc["identifier_exact_duplicate_rows_collapsed"] == 1
    assert qc["identifier_conflicting_duplicate_groups"] == 0
    assert qc["identifier_duplicate_rows_excluded"] == 1


def test_formatter_error_policy_logs_exact_repeat_collapse(
    tmp_path, monkeypatch, capsys,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 100, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "G", "T"],
        "LP": [8.0, 8.0, 4.0],
        "N": [1000.0, 1000.0, 1000.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["magma"],
        duplicate_id_policy="error",
        bcftools=sys.executable,
        overwrite=True,
    ))["magma"]

    assert result["rows_out"] == 2
    assert result["identifier_duplicate_policy"] == "error"
    assert result["identifier_exact_duplicate_rows_collapsed"] == 1
    screen_text = " ".join(capsys.readouterr().out.split())
    assert "Configured duplicate policy : error" in screen_text
    assert "Exact repeated rows collapsed : 1" in screen_text


@pytest.mark.parametrize(
    ("policy", "expected_position", "ranking_column"),
    (
        ("most_significant", 101, "LP"),
        ("highest_maf", 101, "EAF"),
        ("highest_info", 100, "INFO"),
    ),
)
def test_ranked_duplicate_policies_require_one_strict_best_row(
    policy, expected_position, ranking_column,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "LP": [5.0, 9.0, 4.0],
        "EAF": [0.1, 0.4, 0.2],
        "INFO": [0.99, 0.95, 0.90],
    })
    module = load_configuration().modules.formatting

    selected, qc = select_variant_identifiers(
        frame, module, "rsid", duplicate_policy=policy,
    )

    assert selected["POS"].to_list() == [expected_position, 200]
    assert qc["identifier_duplicate_groups_resolved_by_policy"] == 1
    assert qc["identifier_duplicate_groups_unresolved"] == 0
    assert qc["identifier_ranked_winners_retained"] == 1
    assert qc["identifier_duplicate_rows_excluded"] == 1
    assert qc["identifier_duplicate_ranking_column"] == ranking_column


@pytest.mark.parametrize(
    ("policy", "values"),
    (
        ("most_significant", [8.0, 8.0, 4.0]),
        ("highest_maf", [0.2, 0.8, 0.3]),
        ("highest_info", [None, None, 0.9]),
    ),
)
def test_ranked_duplicate_policy_excludes_tied_or_missing_group(policy, values):
    columns = {
        "CHROM": ["1", "1", "2"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "LP": [5.0, 9.0, 4.0],
        "EAF": [0.1, 0.4, 0.3],
        "INFO": [0.95, 0.96, 0.9],
    }
    source = {
        "most_significant": "LP",
        "highest_maf": "EAF",
        "highest_info": "INFO",
    }[policy]
    columns[source] = values
    module = load_configuration().modules.formatting

    selected, qc = select_variant_identifiers(
        pl.DataFrame(columns), module, "rsid", duplicate_policy=policy,
    )

    assert selected["SNP"].to_list() == ["rs2"]
    assert qc["identifier_duplicate_groups_resolved_by_policy"] == 0
    assert qc["identifier_duplicate_groups_unresolved"] == 1
    assert qc["identifier_duplicate_rows_excluded"] == 2


def test_duplicate_policy_cli_overrides_yaml_and_clears_target_overrides(tmp_path):
    from postgwas.modules.formatting.service import _resolved_configuration

    run_config = tmp_path / "formatting.yaml"
    run_config.write_text(
        "modules:\n"
        "  formatting:\n"
        "    variant_identifiers:\n"
        "      default_duplicate_policy: error\n"
        "      target_duplicate_policies: {magma: highest_info}\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args([
        "--vcf", str(FIXTURE),
        "--output-directory", str(tmp_path / "output"),
        "--run-config", str(run_config),
        "--format", "magma", "ldsc",
        "--duplicate-id-policy", "most_significant",
    ])

    configuration = _resolved_configuration(args)

    policy = configuration.modules.formatting.variant_identifiers
    assert policy.default_duplicate_policy == "most_significant"
    assert policy.target_duplicate_policies == {}


def test_default_duplicate_policy_applies_to_every_builtin_format(
    tmp_path, monkeypatch,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "BETA": [0.2, 0.3, -0.1],
        "SE": [0.1, 0.1, 0.2],
        "Z": [2.0, 3.0, -0.5],
        "LP": [5.0, 9.0, 4.0],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
        "N": [1000.0, 1000.0, 1000.0],
        "NEFF": [900.0, 900.0, 900.0],
        "N_CASE": [400.0, 400.0, 400.0],
        "N_CONTROL": [600.0, 600.0, 600.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    for target in TARGETS:
        result = run_formatter_direct(Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(tmp_path / target),
            format=[target],
            bcftools=sys.executable,
            overwrite=True,
        ))[target]

        assert result["identifier_duplicate_policy"] == "exclude_all"
        assert result["identifier_duplicate_rows_excluded"] == 2
        assert result["rows_out"] == 1


def test_multi_format_run_reuses_one_validated_identifier_selection(
    tmp_path, monkeypatch,
):
    from postgwas.modules.formatting import service

    frame = pl.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "BETA": [0.2, 0.3, -0.1],
        "SE": [0.1, 0.1, 0.2],
        "Z": [2.0, 3.0, -0.5],
        "LP": [5.0, 9.0, 4.0],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
        "N": [1000.0, 1000.0, 1000.0],
        "NEFF": [900.0, 900.0, 900.0],
        "N_CASE": [400.0, 400.0, 400.0],
        "N_CONTROL": [600.0, 600.0, 600.0],
    })
    monkeypatch.setattr(service, "load_harmonised_vcf", lambda *args, **kwargs: frame)
    calls = []
    select = service.select_variant_identifiers

    def counted_selection(frame, config, identifier_type, *, duplicate_policy=None):
        calls.append((identifier_type, duplicate_policy))
        return select(
            frame,
            config,
            identifier_type,
            duplicate_policy=duplicate_policy,
        )

    monkeypatch.setattr(service, "select_variant_identifiers", counted_selection)

    results = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=TARGETS,
        bcftools=sys.executable,
        overwrite=True,
    ))

    assert calls == [("rsid", "exclude_all")]
    assert results["magma"]["identifier_selection_reused"] is False
    assert all(
        results[target]["identifier_selection_reused"] is True
        for target in TARGETS[1:]
    )
    assert all(
        result["identifier_uniqueness_validated"] is True
        and result["identifier_validation_rows"] == 1
        and result["identifier_unique_values"] == 1
        for result in results.values()
    )

    written_tables = [
        (results["magma"]["snp_loc_file"], "SNP"),
        (results["magma"]["pval_file"], "SNP"),
        (results["gcta_gene"]["summary_statistics_input_file"], "SNP"),
        (results["susie"]["susie_input"], "SNP"),
        (results["finemap"]["finemap_input"], "rsid"),
        (results["pred_ld"]["files"][1], "snp"),
        (results["ldsc"]["ldsc_file"], "SNP"),
        (results["mixer"]["mixer_input"], "SNP"),
    ]
    for path, identifier in written_tables:
        output = pl.read_csv(path, separator="\t")
        assert output.height == 1
        assert not output[identifier].is_duplicated().any()

    log_text = Path(results["magma"]["log_file"]).read_text(encoding="utf-8")
    assert log_text.count("formatter_unique_identifiers status=PASSED") == len(
        TARGETS
    )
    assert "target=magma" in log_text
    assert "selection_reused=false" in log_text
    assert "target=finemap" in log_text
    assert "selection_reused=true" in log_text


def test_identifier_cache_key_respects_per_target_duplicate_policy(
    tmp_path, monkeypatch,
):
    from postgwas.modules.formatting import service

    run_config = tmp_path / "formatting.yaml"
    run_config.write_text(
        "modules:\n"
        "  formatting:\n"
        "    variant_identifiers:\n"
        "      target_duplicate_policies: {finemap: highest_info}\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "BETA": [0.2, 0.3, -0.1],
        "SE": [0.1, 0.1, 0.2],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
        "LP": [5.0, 9.0, 4.0],
        "N": [1000.0, 1000.0, 1000.0],
        "NEFF": [900.0, 900.0, 900.0],
    })
    monkeypatch.setattr(service, "load_harmonised_vcf", lambda *args, **kwargs: frame)
    calls = []
    select = service.select_variant_identifiers

    def counted_selection(frame, config, identifier_type, *, duplicate_policy=None):
        calls.append((identifier_type, duplicate_policy))
        return select(
            frame,
            config,
            identifier_type,
            duplicate_policy=duplicate_policy,
        )

    monkeypatch.setattr(service, "select_variant_identifiers", counted_selection)

    results = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "output"),
        format=["magma", "finemap"],
        run_config=str(run_config),
        bcftools=sys.executable,
        overwrite=True,
    ))

    assert calls == [
        ("rsid", "exclude_all"),
        ("rsid", "highest_info"),
    ]
    assert results["magma"]["rows_out"] == 1
    assert results["finemap"]["rows_out"] == 2
    assert results["magma"]["identifier_selection_reused"] is False
    assert results["finemap"]["identifier_selection_reused"] is False


def test_uniqueness_invariant_stops_before_formatter_artifacts_are_written(
    tmp_path, monkeypatch,
):
    from postgwas.modules.formatting import service

    frame = pl.DataFrame({
        "CHROM": ["1", "1"],
        "POS": [100, 101],
        "ID": ["rs1", "rs2"],
        "REF": ["A", "C"],
        "ALT": ["G", "T"],
        "LP": [5.0, 4.0],
        "N": [1000.0, 1000.0],
    })
    monkeypatch.setattr(service, "load_harmonised_vcf", lambda *args, **kwargs: frame)
    monkeypatch.setattr(
        service,
        "select_variant_identifiers",
        lambda *args, **kwargs: (
            frame.with_columns(pl.lit("rs1").alias("SNP")),
            {},
        ),
    )

    with pytest.raises(
        FormattingError,
        match=r"identifier validation failed.*2 duplicated rsid.*Refusing to write",
    ):
        run_formatter_direct(Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(tmp_path),
            format=["magma"],
            bcftools=sys.executable,
            overwrite=True,
        ))

    assert not (tmp_path / "STUDY_magma_snp_loc.tsv").exists()
    assert not (tmp_path / "STUDY_magma_p_values.tsv").exists()
    assert "Formatter identifier validation failed" in (
        tmp_path / "logs" / "STUDY_formatter.log"
    ).read_text(encoding="utf-8")


def test_magma_formatter_excludes_all_duplicate_identifier_rows(
    tmp_path, monkeypatch, capsys,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "2", "3", "4"],
        "POS": [100, 500, 200, 300],
        "ID": ["rs1", "rs1", "rs2", "."],
        "REF": ["A", "C", "G", "T"],
        "ALT": ["G", "T", "A", "C"],
        "LP": [5.0, 9.0, 4.0, 3.0],
        "N": [1000.0, 1000.0, 1000.0, 1000.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["magma"],
        bcftools=sys.executable,
        overwrite=True,
    ))["magma"]

    p_values = pl.read_csv(result["pval_file"], separator="\t")
    locations = pl.read_csv(result["snp_loc_file"], separator="\t")
    assert p_values["SNP"].to_list() == ["rs2"]
    assert locations["SNP"].to_list() == ["rs2"]
    assert result["rows_in"] == 4
    assert result["rows_out"] == 1
    assert result["rows_excluded"] == 3
    assert result["identifier_duplicate_groups"] == 1
    assert result["identifier_duplicate_rows"] == 2
    assert result["identifier_duplicate_rows_excluded"] == 2
    assert result["identifier_missing_rows_excluded"] == 1
    log_text = Path(result["log_file"]).read_text(encoding="utf-8")
    assert "WARNING  duplicated rsid identifiers" in log_text
    assert "reason=duplicate_variant_identifier" in log_text
    screen_text = " ".join(capsys.readouterr().out.split())
    assert "Formatter plan · prepare files required by downstream analyses" in screen_text
    assert (
        "Create validated PostGWAS intermediate files for this tool; the "
        "downstream analysis is not run in this step"
    ) in screen_text
    assert "MAGMA input files prepared" in screen_text
    assert "Validate and read the harmonised GWAS-VCF" in screen_text
    assert "Create and validate MAGMA input files" in screen_text
    assert "Total variants in input GWAS-VCF" in screen_text
    assert "MAGMA files created" in screen_text
    assert "Variants read from GWAS-VCF" in screen_text
    assert "Variants written to MAGMA files" in screen_text
    assert "Variants excluded before writing" in screen_text
    assert "Duplicate rsID handling" in screen_text
    assert "1 group containing 2 rows (50.00% of GWAS-VCF rows)" in screen_text
    assert "Resolved by duplicate policy" in screen_text
    assert "0 groups (0%)" in screen_text
    assert "Unresolved groups excluded" in screen_text
    assert "1 group (100%)" in screen_text
    assert "Final identifier validation" in screen_text
    assert "Passed — all 1 retained identifiers are unique" in screen_text
    assert "Rows removed by duplicate policy" in screen_text
    assert "Prepared scientific files" in screen_text
    assert Path(result["snp_loc_file"]).name in screen_text
    assert Path(result["pval_file"]).name in screen_text


def test_pipeline_formatter_keeps_progress_but_suppresses_duplicate_reports(
    tmp_path, monkeypatch, capsys,
):
    frame = pl.DataFrame({
        "CHROM": ["1"],
        "POS": [100],
        "ID": ["rs1"],
        "REF": ["A"],
        "ALT": ["G"],
        "LP": [5.0],
        "N": [1000.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )
    context = {}

    run_formatter_direct(
        Namespace(
            vcf=str(FIXTURE),
            dataset_id="STUDY",
            output_directory=str(tmp_path),
            format=["magma"],
            bcftools=sys.executable,
            overwrite=True,
        ),
        context,
    )

    screen_text = " ".join(capsys.readouterr().out.split())
    assert "Formatter preparation progress" in screen_text
    assert "Total variants in input GWAS-VCF : 1" in screen_text
    assert "Variants written : 1" in screen_text
    assert "Formatter plan · prepare files required by downstream analyses" not in (
        screen_text
    )
    assert "MAGMA input files prepared" not in screen_text
    assert context["formatter"]["magma"]["rows_out"] == 1


def test_magma_formatter_duplicate_policy_is_schema_validated(tmp_path):
    module = load_configuration().modules.formatting
    assert module.variant_identifiers.default_duplicate_policy == "exclude_all"
    assert module.variant_identifiers.target_duplicate_policies == {}

    supported = tmp_path / "supported_duplicate_policies.yaml"
    supported.write_text(
        "variant_identifiers:\n"
        "  target_duplicate_policies:\n"
        "    magma: most_significant\n"
        "    gcta_gene: highest_maf\n"
        "    susie: highest_info\n"
        "    finemap: error\n"
        "    pred_ld: exclude_all\n"
        "    ldsc: most_significant\n"
        "    mixer: highest_info\n",
        encoding="utf-8",
    )
    resolved = load_module_configuration("formatting", supported)
    assert set(resolved.variant_identifiers.target_duplicate_policies) == set(
        TARGETS
    )

    invalid = tmp_path / "invalid_duplicate_policy.yaml"
    invalid.write_text(
        "variant_identifiers:\n"
        "  target_duplicate_policies: {magma: keep_first}\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ConfigurationError,
        match=(
            r"modules\.formatting\.variant_identifiers\."
            r"target_duplicate_policies\.magma"
        ),
    ):
        load_module_configuration("formatting", invalid)


def test_direct_ldsc_without_reference_warns_and_excludes_duplicate_rsids(
    tmp_path, monkeypatch, capsys,
):
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "1"],
        "POS": [100, 101, 200],
        "ID": ["rs1", "rs1", "rs2"],
        "REF": ["A", "A", "C"],
        "ALT": ["G", "C", "T"],
        "Z": [2.0, 3.0, -2.0],
        "LP": [4.0, 6.0, 5.0],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
        "N_CASE": [400.0, 400.0, 400.0],
        "N_CONTROL": [600.0, 600.0, 600.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )

    result = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path),
        format=["ldsc"],
        bcftools=sys.executable,
        overwrite=True,
    ))["ldsc"]

    output = pl.read_csv(result["ldsc_file"], separator="\t")
    assert output["SNP"].to_list() == ["rs2"]
    assert result["rows_in"] == 3
    assert result["rows_out"] == 1
    assert result["rows_excluded"] == 2
    assert result["identifier_duplicate_groups"] == 1
    assert result["identifier_duplicate_rows"] == 2
    assert result["identifier_duplicate_rows_excluded"] == 2
    assert result["identifier_missing_rows_excluded"] == 0
    log_text = Path(result["log_file"]).read_text(encoding="utf-8")
    assert "WARNING  duplicated rsid identifiers" in log_text
    assert "reason=duplicate_variant_identifier" in log_text
    screen_text = " ".join(capsys.readouterr().out.split())
    assert "Rows removed by duplicate policy : 2 (66.67% of GWAS-VCF variants)" in screen_text
    assert "2 rows (66.67% of GWAS-VCF rows)" in screen_text


def test_ldsc_reference_resolves_duplicate_rsid_by_alleles(tmp_path):
    reference = tmp_path / "w_hm3.snplist"
    reference.write_text(
        "SNP A1 A2\nrs1 G A\nrs2 C T\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "SNP": ["rs1", "rs1", "rs2", "rs3"],
        "REF": ["A", "A", "G", "C"],
        "ALT": ["C", "G", "A", "T"],
    })
    module = load_configuration().modules.formatting

    selected, qc = select_ldsc_reference_variants(frame, module, reference)

    assert selected.to_dicts() == [
        {"SNP": "rs1", "REF": "A", "ALT": "G"},
        {"SNP": "rs2", "REF": "G", "ALT": "A"},
    ]
    assert qc == {
        "reference_selection_applied": True,
        "reference_variants": 2,
        "reference_rows_in": 4,
        "reference_rows_out": 2,
        "rows_excluded_not_in_reference": 1,
        "rows_excluded_reference_allele_mismatch": 1,
        "identifier_duplicate_groups": 1,
        "identifier_duplicate_rows": 2,
        "identifier_duplicate_groups_resolved_by_reference": 1,
        "identifier_duplicate_groups_unresolved_after_reference": 0,
    }


def test_ldsc_reference_leaves_ambiguous_group_for_shared_policy(tmp_path):
    reference = tmp_path / "w_hm3.snplist"
    reference.write_text("SNP A1 A2\nrs1 G A\n", encoding="utf-8")
    frame = pl.DataFrame({
        "SNP": ["rs1", "rs1"],
        "REF": ["A", "T"],
        "ALT": ["G", "C"],
    })
    module = load_configuration().modules.formatting

    selected, qc = select_ldsc_reference_variants(frame, module, reference)

    assert selected.height == 2
    assert qc["identifier_duplicate_groups_unresolved_after_reference"] == 1


def test_ldsc_formatter_applies_optional_reference_and_records_provenance(
    tmp_path, monkeypatch, capsys,
):
    from postgwas.modules.formatting import service

    reference = tmp_path / "w_hm3.snplist"
    reference.write_text(
        "SNP A1 A2\nrs1 G A\nrs2 C T\n",
        encoding="utf-8",
    )
    frame = pl.DataFrame({
        "CHROM": ["1", "1", "1", "1"],
        "POS": [100, 101, 200, 300],
        "ID": ["rs1", "rs1", "rs2", "rs3"],
        "REF": ["A", "A", "G", "C"],
        "ALT": ["C", "G", "A", "T"],
        "Z": [2.0, 3.0, -2.0, 1.0],
        "LP": [4.0, 6.0, 5.0, 2.0],
        "EAF": [0.2, 0.3, 0.4, 0.1],
        "INFO": [0.95, 0.96, 0.97, 0.98],
        "N": [1000.0, 1000.0, 1000.0, 1000.0],
        "N_CASE": [400.0, 400.0, 400.0, 400.0],
        "N_CONTROL": [600.0, 600.0, 600.0, 600.0],
    })
    monkeypatch.setattr(service, "load_harmonised_vcf", lambda *args, **kwargs: frame)
    calls = []
    select = service.select_variant_identifiers

    def counted_selection(frame, config, identifier_type, *, duplicate_policy=None):
        calls.append((identifier_type, duplicate_policy))
        return select(
            frame,
            config,
            identifier_type,
            duplicate_policy=duplicate_policy,
        )

    monkeypatch.setattr(service, "select_variant_identifiers", counted_selection)

    results = run_formatter_direct(Namespace(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "output"),
        format=["magma", "ldsc"],
        merge_alleles=str(reference),
        bcftools=sys.executable,
        overwrite=True,
    ))
    result = results["ldsc"]

    assert calls == [("rsid", "exclude_all"), ("rsid", "allow")]
    assert results["magma"]["identifier_selection_reused"] is False
    assert result["identifier_selection_reused"] is False
    assert result["identifier_uniqueness_validated"] is True
    output = pl.read_csv(result["ldsc_file"], separator="\t")
    assert output["SNP"].to_list() == ["rs1", "rs2"]
    assert result["rows_in"] == 4
    assert result["rows_out"] == 2
    assert result["rows_excluded_not_in_reference"] == 1
    assert result["rows_excluded_reference_allele_mismatch"] == 1
    assert result["identifier_duplicate_groups_resolved_by_reference"] == 1
    log_text = Path(result["log_file"]).read_text(encoding="utf-8")
    assert "ldsc_merge_alleles" in log_text
    assert "rows_excluded_not_in_reference=1" in log_text
    raw_screen = capsys.readouterr().out
    screen = " ".join(raw_screen.split())
    assert "LDSC files created" in screen
    assert "Supplied LDSC SNP/allele reference" in screen
    assert "2 unique rsIDs" in screen
    assert "Variants written to LDSC files" in screen
    assert (
        "2 (50.00% of GWAS-VCF variants); 100% of supplied LDSC "
        "SNP/allele reference variants"
    ) in screen
    assert "No rsID match in supplied LDSC SNP/allele reference" in screen
    assert "1 (25.00% of GWAS-VCF variants)" in screen
    assert "Allele mismatch with supplied LDSC SNP/allele reference" in screen
    assert "1 group containing 2 rows (50.00% of GWAS-VCF rows)" in screen
    assert "Resolved by supplied LDSC SNP/allele reference" in screen
    assert "1 group (100% of duplicate groups)" in screen
    assert "Remaining ambiguous groups" in screen
    assert "0 groups (0%)" in screen
    assert "Detailed HTML report" in screen
    duplicate_heading = next(
        line for line in raw_screen.splitlines()
        if "Duplicate rsID handling" in line
    )
    duplicate_count = next(
        line for line in raw_screen.splitlines()
        if "Duplicate groups found" in line
    )
    assert duplicate_heading.startswith("          🔬")
    assert duplicate_count.startswith("              🧮")

    report = Path(result["html_report"])
    html = report.read_text(encoding="utf-8")
    assert "LDSC SNP/allele reference matching" in html
    assert "Supplied LDSC SNP/allele reference variants" in html
    assert "Final LDSC-input coverage of supplied SNP/allele reference" in html
    assert "No rsID match in supplied LDSC SNP/allele reference" in html
    assert "Allele mismatch with supplied LDSC SNP/allele reference" in html
    assert "Groups resolved by supplied LDSC SNP/allele reference" in html
    assert "GWAS-VCF case fraction" in html
    assert "distinct from the chromosome-split LD-score files" in html
    assert "href=\"../STUDY_ldsc_input.tsv\"" in html

    manifest = yaml.safe_load(
        (tmp_path / "output" / "run_metadata" / "formatter_completion.yaml")
        .read_text(encoding="utf-8")
    )
    fingerprint = manifest["input_resources"]["ldsc_merge_alleles"]
    assert Path(fingerprint["path"]) == reference.resolve()
    assert fingerprint["sha256"]


def test_ldsc_formatter_restarts_changed_merge_alleles_then_revalidates(
    tmp_path, monkeypatch,
):
    reference = tmp_path / "w_hm3.snplist"
    reference.write_text("SNP A1 A2\nrs1 G A\n", encoding="utf-8")
    frame = pl.DataFrame({
        "CHROM": ["1"], "POS": [100], "ID": ["rs1"],
        "REF": ["A"], "ALT": ["G"], "Z": [2.0], "LP": [4.0],
        "EAF": [0.2], "INFO": [0.95],
        "N_CASE": [400.0], "N_CONTROL": [600.0],
    })
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.load_harmonised_vcf",
        lambda *args, **kwargs: frame,
    )
    common = dict(
        vcf=str(FIXTURE),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "output"),
        format=["ldsc"],
        merge_alleles=str(reference),
        bcftools=sys.executable,
    )
    run_formatter_direct(Namespace(**common, overwrite=True))
    reference.write_text("SNP A1 A2\nrs1 C A\n", encoding="utf-8")

    with pytest.raises(FormattingError, match="No VCF records match"):
        run_formatter_direct(Namespace(**common))

    log_text = (
        tmp_path / "output" / "logs" / "STUDY_formatter.log"
    ).read_text(encoding="utf-8")
    assert "input resource ldsc_merge_alleles changed" in log_text
    assert "reason=changed_inputs" in log_text


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


def test_ldsc_sample_prevalence_aggregation_is_schema_validated(tmp_path):
    config = load_configuration().modules.formatting
    assert config.ldsc_sample_prevalence.aggregation == "median"

    mean_config = tmp_path / "mean_prevalence_aggregation.yaml"
    mean_config.write_text(
        "ldsc_sample_prevalence:\n  aggregation: mean\n",
        encoding="utf-8",
    )
    assert load_module_configuration(
        "formatting", mean_config,
    ).ldsc_sample_prevalence.aggregation == "mean"

    invalid = tmp_path / "invalid_prevalence_aggregation.yaml"
    invalid.write_text(
        "ldsc_sample_prevalence:\n  aggregation: sum\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ConfigurationError,
        match=r"modules\.formatting\.ldsc_sample_prevalence\.aggregation",
    ):
        load_module_configuration("formatting", invalid)


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
    from postgwas.pipeline import runners

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
    monkeypatch.setattr(runners, "_validate_current_pipeline_vcf", lambda *_: None)
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

    configuration = load_configuration(cli_overrides={
        "modules.magma.input.ld_reference_prefix": str(reference),
    })
    context = _magma_pipeline_context(configuration, reference)
    result = runners.run_formatter_runner(args, context)

    assert result == {"magma": {"ok": True}}
    assert captured["target_types"] == {"magma": "unique"}
    assert "magma" in captured["formats"]
    assert captured["observation"]["unique_ids"] == 2
    assert args.output_directory == str(output)


@pytest.mark.parametrize(
    ("progress_kind", "analysis_module"),
    (
        ("magma", "magma"),
        ("magmacovar", "magmacovar"),
        ("pops", "pops"),
    ),
)
def test_magma_pipeline_preflight_reports_detailed_validated_evidence(
    tmp_path, monkeypatch, progress_kind, analysis_module,
):
    from postgwas.pipeline import runners

    class CapturingProgress:
        def __init__(self):
            self.started = []
            self.completed = {}

        def start(self, number):
            self.started.append(number)

        def complete(self, number, **values):
            self.completed[number] = values

    vcf = tmp_path / "study.vcf.gz"
    vcf.write_text("vcf\n", encoding="utf-8")
    gene_locations = tmp_path / "genes.loc"
    gene_locations.write_text("locations\n", encoding="utf-8")
    pathways = tmp_path / "pathways.gmt"
    pathways.write_text("pathways\n", encoding="utf-8")
    reference = tmp_path / "reference"
    output = tmp_path / "output"
    progress = CapturingProgress()

    monkeypatch.setattr(
        runners,
        "_validate_current_pipeline_vcf",
        lambda *_: {
            "indexed": Namespace(variant_count=100),
            "harmonised": {
                "genome_build": "GRCh37",
                "postgwas_dataset_id": "VCF_STUDY",
            },
        },
    )

    def configure_identifiers(args, *_):
        args.variant_id_types = {"magma": "unique"}
        args.variant_id_observations = {
            "magma": {
                "variant_id_type": "unique",
                "variants": 250,
                "consumers": ["MAGMA"],
            },
        }

    monkeypatch.setattr(
        runners, "configure_reference_variant_identifiers", configure_identifiers,
    )
    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.validate_magma_ld_reference_and_runtime",
        lambda *_: pytest.fail("MAGMA resource validation was repeated"),
    )
    monkeypatch.setattr(
        "postgwas.modules.magma.annotations.read_gene_locations",
        lambda *_: pytest.fail("MAGMA gene-location validation was repeated"),
    )
    monkeypatch.setattr(
        "postgwas.modules.magma.analysis.validate_magma_gene_set_references",
        lambda *_, **__: pytest.fail("MAGMA pathway validation was repeated"),
    )
    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        lambda *_: {"magma": {"ok": True}},
    )
    args = Namespace(
        modules=[analysis_module],
        format=None,
        vcf=str(vcf),
        magma_ld_reference=str(reference),
        gene_location_file=str(gene_locations),
        gene_set_file=str(pathways),
        magma=sys.executable,
        bcftools=sys.executable,
        dataset_id="STUDY",
        output_directory=str(output),
        run_config=None,
        _step_num="01",
        _pipeline_stage_progress=progress,
    )
    if progress_kind == "magmacovar":
        from postgwas.modules.magmacovar.stages import (
            magmacovar_pipeline_progress_plan,
        )

        args._pipeline_progress_plan = magmacovar_pipeline_progress_plan(args)
    elif progress_kind == "pops":
        from postgwas.modules.pops.stages import pops_pipeline_progress_plan

        args._pipeline_progress_plan = pops_pipeline_progress_plan(args)
    else:
        from postgwas.modules.magma.reporting import magma_pipeline_progress_plan

        args._pipeline_progress_plan = magma_pipeline_progress_plan(args)

    gene_set_plans = {
        "positional": {
            "status": "ready",
            "pathway_file": str(pathways),
            "gene_reference_file": str(gene_locations),
            "minimum_overlap": 0.5,
            "input_metadata": {"detected_format": "gmt", "gene_sets": 4},
            "validation": {
                "identifier_resolution": {
                    "input_unique_ids": 10,
                    "location_primary_unique_ids": 3,
                    "location_alternate_unique_ids": 2,
                    "location_ambiguous_alternate_ids": 1,
                    "alternate_identifier_column": 6,
                    "direct_primary_matches": 0,
                    "alternate_candidate_matches": 2,
                    "unmatched_input_ids": 8,
                    "identifier_source": "alternate_gene_id",
                },
            },
        },
    }
    configuration = load_configuration(cli_overrides={
        "run.dataset_id": "STUDY",
        "run.output_directory": str(output),
        "modules.magma.input.ld_reference_prefix": str(reference),
        "modules.magma.input.gene_location_file": str(gene_locations),
        "modules.magma.input.gene_set_file": str(pathways),
    })
    context = _magma_pipeline_context(
        configuration,
        reference,
        include_gene_sets=progress_kind == "magma",
        gene_set_plans=gene_set_plans,
        primary_gene_location_metrics={
            "path": str(gene_locations),
            "genes": 3,
            "alternate_unique_ids": 2,
            "ambiguous_alternate_ids": 1,
        },
    )
    assert runners.run_formatter_runner(args, context) == {"magma": {"ok": True}}
    assert progress.started == (
        [1, 2, 3, 4, 5] if progress_kind == "magma" else [1, 2, 3]
    )
    fields = {
        number: {field[1]: field[2] for field in values["outcome_fields"] if len(field) == 3}
        for number, values in progress.completed.items()
    }
    assert fields[1]["Total variants"] == 100
    assert fields[1]["Genome build inferred from VCF"] == "GRCh37"
    assert fields[1]["VCF embedded dataset/sample"] == "VCF_STUDY"
    assert fields[2]["Total BIM variants"] == 250
    assert fields[2]["Variant identifiers"] == "chromosome-position-allele ID"
    assert "not independently verifiable" in fields[2][
        "Build and population provenance"
    ]
    assert fields[3]["Reference genes"] == 3
    assert fields[3]["Duplicated alternate IDs"] == 1
    if progress_kind != "magma":
        assert set(fields) == {1, 2, 3}
        return
    assert fields[4]["Total pathway records"] == 4
    assert fields[4]["Unique gene identifiers"] == 10
    assert fields[5][
        "Pathway-file IDs found in gene-location column 1"
    ] == "0/10 (0.00%)"
    assert fields[5][
        "Pathway-file IDs found in gene-location column 6"
    ] == "2/10 (20.00%)"
    assert fields[5][
        "Gene-location column 6 IDs found in pathway file"
    ] == "2/2 (100.00%)"
    assert fields[5]["Original pathway identifiers"] == "unchanged"
    assert fields[5]["Pathway analysis"] == "available"


def test_magma_pipeline_analysis_uses_exact_formatter_outputs(
    tmp_path, monkeypatch,
):
    from postgwas.pipeline.runners import run_magma_runner

    formatter_locations = tmp_path / "formatter" / "STUDY_magma_snp_loc.tsv"
    formatter_p_values = tmp_path / "formatter" / "STUDY_magma_p_values.tsv"
    captured = {}

    def fake_magma(args, ctx, **kwargs):
        captured["snp_location_file"] = args.snp_location_file
        captured["p_value_file"] = args.p_value_file
        captured["pipeline_resources"] = kwargs["pipeline_resources"]
        return {"magma_genes_out": "genes.out"}

    monkeypatch.setattr(
        "postgwas.modules.magma.service.run_magma_direct",
        fake_magma,
    )
    args = Namespace(
        output_directory=str(tmp_path / "pipeline"),
        _step_num="02",
    )
    context = {
        "formatter": {
            "magma": {
                "snp_loc_file": str(formatter_locations),
                "pval_file": str(formatter_p_values),
            },
        },
    }

    result = run_magma_runner(args, context)

    assert result == {"magma_genes_out": "genes.out"}
    assert captured == {
        "snp_location_file": str(formatter_locations),
        "p_value_file": str(formatter_p_values),
        "pipeline_resources": None,
    }
    assert context["magma"] == result
    assert args.output_directory == str(tmp_path / "pipeline")


def test_heritability_pipeline_passes_merge_alleles_to_ldsc_formatter(
    tmp_path, monkeypatch,
):
    from postgwas.pipeline.runners import run_formatter_runner

    reference = tmp_path / "w_hm3.snplist"
    reference.write_text("SNP A1 A2\nrs1 G A\n", encoding="utf-8")
    captured = {}

    def fake_formatter(args, ctx):
        captured["formats"] = list(args.format)
        captured["target_types"] = dict(args.variant_id_types)
        captured["merge_alleles"] = args.merge_alleles
        return {"ldsc": {"ldsc_file": "formatted.tsv"}}

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        fake_formatter,
    )
    output = tmp_path / "output"
    args = Namespace(
        modules=["heritability"],
        vcf=str(FIXTURE),
        merge_alleles=str(reference),
        output_directory=str(output),
        dataset_id="STUDY",
        bcftools=sys.executable,
        _step_num="01",
    )

    result = run_formatter_runner(args, {})

    assert result == {"ldsc": {"ldsc_file": "formatted.tsv"}}
    assert captured == {
        "formats": ["ldsc"],
        "target_types": {"ldsc": "rsid"},
        "merge_alleles": str(reference),
    }
    assert args.output_directory == str(output)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is required")
def test_pipeline_formatter_restarts_a_manifest_with_no_declared_outputs(tmp_path):
    from postgwas.pipeline.runners import run_formatter_runner

    output = tmp_path / "pipeline"
    args = Namespace(
        modules=["formatter"],
        vcf=str(FIXTURE),
        format=["gcta_gene"],
        output_directory=str(output),
        dataset_id="STUDY",
        bcftools=shutil.which("bcftools"),
        resume=True,
        _step_num="01",
    )
    first = run_formatter_runner(args, {})
    artifact = Path(first["gcta_gene"]["summary_statistics_input_file"])
    artifact.unlink()

    restarted = run_formatter_runner(args, {})

    assert artifact.is_file()
    assert restarted["gcta_gene"].get("resumed") is not True
    assert args.output_directory == str(output)
    log = output / "01_formatter" / "logs" / "STUDY_formatter.log"
    assert "reason=incomplete_outputs" in log.read_text(
        encoding="utf-8"
    )


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
    assert result["sample_prevalence_aggregation"] is None
    assert result["sample_prevalence_variants"] == 0
    assert result["sample_prevalence_minimum"] is None
    assert result["sample_prevalence_maximum"] is None
    assert result["sample_prevalence_case_count_minimum"] is None
    assert result["sample_prevalence_case_count_maximum"] is None
    assert result["sample_prevalence_control_count_minimum"] is None
    assert result["sample_prevalence_control_count_maximum"] is None


def test_ldsc_sample_prevalence_supports_configured_median_and_mean(tmp_path):
    frame = pl.DataFrame({
        "SNP": ["rs1", "rs2", "rs3"],
        "ALT": ["G", "T", "A"],
        "REF": ["A", "C", "G"],
        "Z": [2.0, -2.0, 3.0],
        "LP": [8.0, 6.0, 10.0],
        "N_CASE": [100.0, 400.0, 900.0],
        "N_CONTROL": [900.0, 600.0, 100.0],
        "EAF": [0.2, 0.3, 0.4],
        "INFO": [0.95, 0.96, 0.97],
    })
    default_config = load_configuration().modules.formatting

    median_result = export_ldsc(
        frame, tmp_path / "median", "STUDY", default_config, overwrite=False,
    )
    mean_yaml = tmp_path / "mean.yaml"
    mean_yaml.write_text(
        "ldsc_sample_prevalence:\n  aggregation: mean\n",
        encoding="utf-8",
    )
    mean_config = load_module_configuration("formatting", mean_yaml)
    mean_result = export_ldsc(
        frame, tmp_path / "mean", "STUDY", mean_config, overwrite=False,
    )

    assert median_result["sample_prev"] == pytest.approx(0.4)
    assert median_result["sample_prevalence_aggregation"] == "median"
    assert mean_result["sample_prev"] == pytest.approx((0.1 + 0.4 + 0.9) / 3)
    assert mean_result["sample_prevalence_aggregation"] == "mean"
    for result in (median_result, mean_result):
        assert result["sample_prevalence_variants"] == 3
        assert result["sample_prevalence_minimum"] == pytest.approx(0.1)
        assert result["sample_prevalence_maximum"] == pytest.approx(0.9)
        assert result["sample_prevalence_case_count_minimum"] == 100
        assert result["sample_prevalence_case_count_maximum"] == 900
        assert result["sample_prevalence_control_count_minimum"] == 100
        assert result["sample_prevalence_control_count_maximum"] == 900


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
