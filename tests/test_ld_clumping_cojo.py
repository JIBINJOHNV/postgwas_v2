from argparse import Namespace
import csv
from pathlib import Path
import shutil

import polars as pl
import pytest

from postgwas.config import load_configuration
from postgwas.core.contracts import Artifact, ModuleResult, RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.modules.ld_clumping.ld_prune_cojo import (
    group_cojo_selected_signals,
    prepare_cojo_exclusion_list,
)
from postgwas.modules.ld_clumping.service import (
    LDClumpingInputValidation,
    LDClumpingPreflight,
    _owned_analysis_paths,
    _validate_cojo_reference,
    run_ld_clump_direct,
)
from postgwas.modules.gcta_cojo.service import GctaCojoReferenceValidation
from postgwas.pipeline.runners import run_formatter_runner


def _cojo_frame(rows):
    return pl.DataFrame(
        rows,
        schema=[
            "Chr", "SNP", "bp", "freq", "refA", "b", "se", "p", "n",
            "freq_geno", "bJ", "bJ_se", "pJ", "LD_r",
            "estimation_status",
        ],
        orient="row",
    )


def _row(chromosome, snp, position, p_value, joint_p):
    return (
        chromosome, snp, position, 0.2, "A", 0.2, 0.05, p_value, 10000,
        0.21, 0.15, 0.05, joint_p, 0.1, "estimated",
    )


def test_cojo_grouping_uses_consecutive_250kb_gap_and_preserves_all_signals(
    tmp_path,
):
    configuration = load_configuration()
    ld_module = configuration.modules.ld_clumping.model_copy(deep=True)
    ld_module.remove_mhc = False
    cojo_module = configuration.modules.gcta_cojo
    source = tmp_path / "normalized.tsv"
    _cojo_frame([
        _row("1", "rs1", 100, 1e-9, 1e-8),
        _row("1", "rs2", 250100, 1e-10, 1e-7),
        _row("1", "rs3", 500101, 1e-8, 1e-9),
        _row("2", "rs4", 100, 1e-12, 1e-11),
    ]).write_csv(source, separator="\t")

    result = group_cojo_selected_signals(
        normalized_result=source,
        selected_signals_path=tmp_path / "signals.tsv",
        loci_path=tmp_path / "loci.tsv",
        ld_module=ld_module,
        cojo_module=cojo_module,
    )

    assert result["selected_signals"] == 4
    assert result["genomic_loci"] == 3
    signals = pl.read_csv(result["selected_signals_file"], separator="\t")
    assert signals["COJO_locus"].to_list() == [1, 1, 2, 3]
    assert signals["is_locus_index"].to_list() == [False, True, True, True]
    loci = pl.read_csv(result["loci_file"], separator="\t")
    assert loci["index_SNP"].to_list() == ["rs2", "rs3", "rs4"]


def test_cojo_grouping_writes_valid_empty_tables(tmp_path):
    configuration = load_configuration()
    ld_module = configuration.modules.ld_clumping.model_copy(deep=True)
    ld_module.remove_mhc = False
    cojo_module = configuration.modules.gcta_cojo
    source = tmp_path / "normalized.tsv"
    _cojo_frame([]).write_csv(source, separator="\t")

    result = group_cojo_selected_signals(
        normalized_result=source,
        selected_signals_path=tmp_path / "signals.tsv",
        loci_path=tmp_path / "loci.tsv",
        ld_module=ld_module,
        cojo_module=cojo_module,
    )

    assert result["selected_signals"] == 0
    assert result["genomic_loci"] == 0
    assert pl.read_csv(result["selected_signals_file"], separator="\t").is_empty()
    assert pl.read_csv(result["loci_file"], separator="\t").is_empty()


def test_cojo_genomic_control_preserves_gcta_pvalue_labels(tmp_path):
    configuration = load_configuration()
    ld_module = configuration.modules.ld_clumping.model_copy(deep=True)
    ld_module.remove_mhc = False
    cojo_module = configuration.modules.gcta_cojo.model_copy(deep=True)
    cojo_module.analysis = cojo_module.analysis.model_copy(update={
        "genomic_control": True,
    })
    source = tmp_path / "normalized.tsv"
    _cojo_frame([
        _row("1", "rs1", 100, 1e-9, 1e-8),
    ]).rename({"p": "p_GC", "pJ": "pJ_GC"}).write_csv(
        source, separator="\t",
    )

    result = group_cojo_selected_signals(
        normalized_result=source,
        selected_signals_path=tmp_path / "signals.tsv",
        loci_path=tmp_path / "loci.tsv",
        ld_module=ld_module,
        cojo_module=cojo_module,
    )

    signal_table = result["_report_tables"]["cojo_selected_signals"]
    locus_table = result["_report_tables"]["cojo_loci"]
    assert "p_GC" in signal_table["columns"]
    assert "pJ_GC" in signal_table["columns"]
    assert "p" not in signal_table["columns"]
    assert "pJ" not in signal_table["columns"]
    assert "index_p_GC" in locus_table["columns"]
    assert "index_pJ_GC" in locus_table["columns"]
    locus_columns = pl.read_csv(result["loci_file"], separator="\t").columns
    assert "index_p_GC" in locus_columns
    assert "index_pJ_GC" in locus_columns


def test_cojo_exclusion_list_merges_configured_ids_with_reference_mhc(tmp_path):
    configuration = load_configuration()
    ld_module = configuration.modules.ld_clumping
    cojo_module = configuration.modules.gcta_cojo.model_copy(deep=True)
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bim").write_text(
        "6 rsMHC 0 30000000 A G\n"
        "6 rsOutside 0 40000000 C T\n",
        encoding="utf-8",
    )
    configured = tmp_path / "configured.txt"
    configured.write_text("rsCustom\nrsMHC\n", encoding="utf-8")
    cojo_module.inputs = cojo_module.inputs.model_copy(update={
        "exclude_snps": str(configured),
    })

    result = prepare_cojo_exclusion_list(
        reference_prefix=prefix,
        destination=tmp_path / "combined.txt",
        ld_module=ld_module,
        cojo_module=cojo_module,
    )

    assert result["configured_exclusions"] == 2
    assert result["mhc_exclusions"] == 1
    assert result["combined_exclusions"] == 2
    assert Path(result["path"]).read_text(encoding="utf-8").splitlines() == [
        "rsCustom", "rsMHC",
    ]


def test_cojo_exclusion_list_accepts_reference_without_mhc_variants(tmp_path):
    configuration = load_configuration()
    ld_module = configuration.modules.ld_clumping
    cojo_module = configuration.modules.gcta_cojo
    prefix = tmp_path / "chr1_reference"
    Path(str(prefix) + ".bim").write_text(
        "1 rs1 0 100 A G\n",
        encoding="utf-8",
    )

    result = prepare_cojo_exclusion_list(
        reference_prefix=prefix,
        destination=tmp_path / "combined.txt",
        ld_module=ld_module,
        cojo_module=cojo_module,
    )

    assert result["path"] is None
    assert result["mhc_exclusions"] == 0
    assert result["combined_exclusions"] == 0


def test_cojo_reference_preflight_forwards_canonical_logger(
    monkeypatch, tmp_path,
):
    configuration = load_configuration(cli_overrides={
        "modules.gcta_cojo.mode": "slct",
    })
    reference = type("Reference", (), {
        "prefix": tmp_path / "reference",
        "files": {"bim": tmp_path / "reference.bim"},
    })()
    observed = {}
    logger = object()

    def validate(configuration, *, logger=None):
        observed["logger"] = logger
        return reference

    monkeypatch.setattr(
        "postgwas.modules.gcta_cojo.service.validate_gcta_cojo_reference",
        validate,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service.configure_reference_variant_identifiers",
        lambda *args, **kwargs: None,
    )

    result = _validate_cojo_reference(
        Namespace(),
        configuration,
        logger=logger,
    )

    assert result is reference
    assert observed["logger"] is logger


def test_cojo_workspace_uses_flat_shared_hierarchy_and_scoped_filenames():
    module = load_configuration().modules.ld_clumping
    layout = module.output_layout

    assert layout.cojo_formatter_directory == "cojo/formatter"
    assert layout.cojo_root_directory == "cojo"
    assert module.cojo.subordinate_dataset_id.format(
        dataset_id="STUDY", population="EUR",
    ) == "STUDY_EUR"


def test_shared_cojo_directories_are_not_outer_quarantine_targets(tmp_path):
    module = load_configuration().modules.ld_clumping

    paths = _owned_analysis_paths(tmp_path, module, "STUDY")

    assert tmp_path / "cojo" not in paths
    assert tmp_path / "cojo" / "formatter" not in paths


@pytest.mark.parametrize(
    ("key", "value", "message"),
    (
        (
            "modules.ld_clumping.cojo.subordinate_dataset_id",
            "{dataset_id}",
            "must contain {dataset_id} and {population} exactly once",
        ),
        (
            "modules.ld_clumping.cojo.formatter_resolved_config_file",
            "run_metadata/resolved_config.yaml",
            "must contain {dataset_id} exactly once",
        ),
        (
            "modules.ld_clumping.cojo.gcta_resolved_config_file",
            "run_metadata/{dataset_id}_resolved_config.yaml",
            "must contain {dataset_id} and {mode} exactly once",
        ),
    ),
)
def test_cojo_shared_hierarchy_rejects_unscoped_subordinate_names(
    key, value, message,
):
    with pytest.raises(ConfigurationError, match=message):
        load_configuration(cli_overrides={key: value})


def test_ld_service_reports_cojo_signals_and_loci_in_html(
    monkeypatch, tmp_path, capsys,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf")
    output = tmp_path / "output"
    prefix = tmp_path / "reference"
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.methods": ["cojo-slct"],
        "modules.ld_clumping.inputs.vcf": str(vcf),
        "modules.ld_clumping.inputs.dataset_id": "STUDY",
        "modules.ld_clumping.output_directory": str(output),
        "modules.ld_clumping.remove_mhc": False,
        "modules.gcta_cojo.genome_build": "GRCh37",
        "modules.gcta_cojo.reference.prefix": str(prefix),
        "modules.gcta_cojo.reference.population": "EUR",
        "run.dataset_id": "STUDY",
        "run.output_directory": str(output),
        "logging.show_screen": True,
        "logging.show_progress": True,
    })
    preflight = LDClumpingPreflight(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        tabix=None,
        reference_directory=None,
        reference_manifest=None,
        gcta="gcta64",
        cojo_reference_prefix=prefix,
        cojo_reference_files={
            suffix: Path(str(prefix) + "." + suffix)
            for suffix in ("bed", "bim", "fam")
        },
        cojo_reference_validation=GctaCojoReferenceValidation(
            prefix=prefix,
            files={
                suffix: Path(str(prefix) + "." + suffix)
                for suffix in ("bed", "bim", "fam")
            },
            samples=503,
            executable="gcta64",
            version="1.94.1",
            software={"path": "gcta64"},
        ),
        cojo_reference_identifier_observation={
            "variant_id_type": "unique",
            "variants": 13_457_498,
        },
    )
    input_validation = LDClumpingInputValidation(
        configuration=configuration,
        vcf=vcf,
        output_directory=output,
        dataset_id="STUDY",
        bcftools="bcftools",
        genome_build="GRCh37",
        sample="STUDY",
        variant_count=100,
        contigs=("1",),
        required_fields=("FORMAT/LP",),
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service._validate_ld_clumping_input",
        lambda *args, **kwargs: input_validation,
    )
    monkeypatch.setattr(
        "postgwas.modules.ld_clumping.service._validate_ld_clumping_references",
        lambda *args, **kwargs: preflight,
    )
    subordinate = {}

    def formatter(*args, configuration, **kwargs):
        module = configuration.modules.formatting
        subordinate["formatter_dataset_id"] = configuration.run.dataset_id
        subordinate["formatter_resolved_config"] = (
            module.runtime.resolved_config_file
        )
        subordinate["formatter_completion"] = (
            module.runtime.completion_manifest_file
        )
        output_pattern = module.exports["gcta_gene"].outputs[
            "summary_statistics"
        ].output_file
        path = Path(configuration.run.output_directory) / output_pattern.format(
            dataset_id=configuration.run.dataset_id,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "SNP A1 A2 freq b se p N\nrs1 A G 0.2 0.2 0.05 1e-9 10000\n",
            encoding="utf-8",
        )
        return {"gcta_gene": {"summary_statistics_input_file": str(path)}}

    def gcta(*args, configuration, **kwargs):
        module = configuration.modules.gcta_cojo
        subordinate["gcta_output_contract"] = kwargs[
            "parallel_slct_output_contract"
        ]
        subordinate["gcta_dataset_id"] = configuration.run.dataset_id
        subordinate["gcta_resolved_config"] = (
            module.output_layout.resolved_config_file
        )
        values = {
            "dataset_id": configuration.run.dataset_id,
            "mode": module.mode,
        }
        normalized = Path(configuration.run.output_directory) / (
            module.output_layout.normalized_result.format(**values)
        )
        log = Path(configuration.run.output_directory) / (
            module.output_layout.log_file.format(**values)
        )
        normalized.parent.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        _cojo_frame([
            _row("1", "rs1", 100, 1e-9, 1e-8),
            _row("1", "rs2", 200000, 1e-10, 1e-9),
        ]).write_csv(normalized, separator="\t")
        log.write_text("GCTA completed\n", encoding="utf-8")
        return ModuleResult(
            "gcta_cojo",
            artifacts={
                "normalized_results": Artifact("gcta", normalized),
                "gcta_log": Artifact("log", log),
            },
            metrics={
                "summary_variants": 100,
                "reference_overlap_variants": 90,
                "reference_overlap_fraction": 0.9,
                "reference_missing_variants": 10,
                "reference_allele_mismatch_variants": 0,
                "reference_samples": 5000,
                "gcta_version": "1.94.1",
                "execution_strategy": "chromosome_parallel",
                "parallel_output_contract": "selection_only",
                "conditional_reconstruction_status": "not_requested",
                "command": ["gcta64", "--cojo-slct", "--cojo-wind", "10000"],
            },
            warnings=(
                "The LD reference has 5000 samples; review panel suitability.",
                "10 of 100 GWAS variants are absent from the LD-reference BIM "
                "and will not enter COJO.",
            ),
        )

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct", formatter,
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_cojo.service.run_gcta_cojo_direct", gcta,
    )

    result = run_ld_clump_direct(
        type("Args", (), {})(), configuration=configuration,
    )

    screen = capsys.readouterr().out
    compact_screen = " ".join(screen.split())
    assert result["ld_clump_cojo"]["selected_signals"] == 2
    assert result["ld_clump_cojo"]["genomic_loci"] == 1
    cojo_root = output / "cojo"
    assert result["ld_clump_cojo"]["cojo_output_directory"] == str(cojo_root)
    assert result["ld_clump_cojo"]["formatter_output"] == str(
        cojo_root / "formatter" / "STUDY_EUR_gcta.ma"
    )
    assert result["ld_clump_cojo"]["output_files"][
        "gcta_normalized_results"
    ] == str(cojo_root / "results" / "STUDY_EUR_slct_gcta_cojo.tsv")
    assert subordinate == {
        "formatter_dataset_id": "STUDY_EUR",
        "formatter_resolved_config": (
            "run_metadata/{dataset_id}_formatter_resolved_config.yaml"
        ),
        "formatter_completion": (
            "run_metadata/{dataset_id}_formatter_completion.yaml"
        ),
        "gcta_dataset_id": "STUDY_EUR",
        "gcta_resolved_config": (
            "run_metadata/{dataset_id}_{mode}_gcta_cojo_resolved_config.yaml"
        ),
        "gcta_output_contract": "selection_only",
    }
    assert "GCTA-COJO stepwise method" in screen
    assert "Completed 1/4 · Validate the input GWAS-VCF" in screen
    assert "Input summary-statistics VCF" in screen
    assert "Total variants" in screen
    assert "Completed 2/4 · Validate selected LD-reference inputs" in screen
    assert "PLINK LD reference for GCTA-COJO" in screen
    assert "BED, BIM, FAM present and non-empty" in compact_screen
    assert "Total BIM variants" in screen
    assert "13,457,498" in screen
    assert "LD-reference samples" in screen
    assert "503" in screen
    assert "chromosome" in screen
    assert "allele ID" in compact_screen
    assert "BIM structural validation" in screen
    assert "Build and population provenance" in screen
    assert "not independently verifiable" in compact_screen
    assert "COJO execution plan" not in screen
    assert "Selection strategy" not in screen
    assert "Available compute" not in screen
    assert "Scheduling rule" not in screen
    assert "GWAS variants matched to LD reference" in screen
    assert "90 of 100 · 90.00% matched by SNP ID and allele pair" in screen
    assert "COJO warnings" in screen
    assert "Warning 1" in screen
    assert "review panel suitability" in screen
    assert "Warning 2" in screen
    assert "will not enter COJO" in compact_screen
    assert "COJO-selected signals" in screen
    assert "Physical loci" in screen
    assert "COJO output scope" in screen
    assert "full-genome conditional table (.cma) not generated" in compact_screen
    canonical_log = Path(result["canonical_log"]).read_text(encoding="utf-8")
    assert "GCTA-COJO stepwise method" in canonical_log
    assert "cojo_execution_setup" in canonical_log
    assert "--cojo-slct" in canonical_log
    with Path(result["summary_csv"]).open(encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
    warning_rows = [
        row for row in summary_rows if row["record_type"] == "cojo_warning"
    ]
    assert len(warning_rows) == 2
    assert warning_rows[0]["detail_name"] == "warning_1"
    assert "review panel suitability" in warning_rows[0]["detail_value"]
    assert warning_rows[1]["detail_name"] == "warning_2"
    assert "will not enter COJO" in warning_rows[1]["detail_value"]
    assert warning_rows[0]["cojo_reference_overlap_fraction"] == "0.9"
    method_row = next(
        row for row in summary_rows
        if row["record_type"] == "method" and row["method"] == "cojo-slct"
    )
    assert method_row["cojo_parallel_output_contract"] == "selection_only"
    assert method_row["cojo_configured_parallel_output_contract"] == (
        "selection_only"
    )
    assert method_row["cojo_conditional_reconstruction_status"] == (
        "not_requested"
    )
    html = Path(result["html_report"]).read_text(encoding="utf-8")
    assert "COJO-selected signals" in html
    assert "rs1" in html
    assert "rs2" in html
    assert "GCTA-COJO model and physical locus policy" in html
    assert "GWAS variants matched to LD reference" in html
    assert "90 of 100 · 90.00% matched by SNP ID and allele pair" in html
    assert "COJO completed with 2 warnings" in html
    assert "COJO warnings" in html
    assert "Reason and consequence" in html
    assert "review panel suitability" in html
    assert "will not enter COJO" in html
    assert "--cojo-wind" in html
    assert "COJO output scope" in html
    assert "full-genome" in html
    assert "--cojo-cond" in html
    assert "not used to select signals or define the reported physical loci" in html


def test_ld_pipeline_formatter_adds_gcta_input_for_cojo_method(
    monkeypatch, tmp_path,
):
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bim").write_text(
        "1 1_100_A_G 0 100 G A\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    captured = {}

    def formatter(args, ctx):
        captured["formats"] = list(args.format)
        captured["identifier_types"] = dict(args.variant_id_types)
        return {"gcta_gene": {"summary_statistics_input_file": "study.ma"}}

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct", formatter,
    )
    monkeypatch.setattr(
        "postgwas.pipeline.runners._validate_current_pipeline_vcf",
        lambda *_args: None,
    )
    output = tmp_path / "output"
    args = Namespace(
        modules=["ld_clump"],
        clumping_methods=["cojo-slct"],
        cojo_reference_prefix=str(prefix),
        vcf=str(vcf),
        output_directory=str(output),
        dataset_id="STUDY",
        bcftools=shutil.which("bcftools") or shutil.which("python"),
        _step_num="01",
    )
    configuration = load_configuration(cli_overrides={
        "modules.ld_clumping.methods": ["cojo-slct"],
        "modules.gcta_cojo.reference.prefix": str(prefix),
        "modules.gcta_cojo.reference.population": "EUR",
        "modules.gcta_cojo.genome_build": "GRCh37",
    })
    reference = GctaCojoReferenceValidation(
        prefix=prefix,
        files={"bim": Path(str(prefix) + ".bim")},
        samples=1,
        executable="gcta64",
        version="1.94.1",
        software={"path": "gcta64"},
    )
    ld_resources = Namespace(
        configuration=configuration,
        cojo_reference_validation=reference,
    )
    context = RunContext(validations={
        "ld_clump": PipelinePreflightEvidence(
            module="ld_clump",
            input_vcf={},
            resources=ld_resources,
        ),
    })

    result = run_formatter_runner(args, context)

    assert result["gcta_gene"]["summary_statistics_input_file"] == "study.ma"
    assert captured["formats"] == ["gcta_gene"]
    assert captured["identifier_types"] == {"gcta_gene": "unique"}
    assert args.output_directory == str(output)
