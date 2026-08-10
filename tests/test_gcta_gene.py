import argparse
from argparse import Namespace
from pathlib import Path

import pytest
import polars as pl
import yaml

from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.ui import format_cli_default
from postgwas.modules.gcta_gene import cli as gcta_gene_cli
from postgwas.modules.gcta_gene.adapters import build_gcta_command
from postgwas.modules.gcta_gene.cli import build_parser
from postgwas.modules.gcta_gene.errors import GctaGeneError
from postgwas.modules.gcta_gene.reporting import build_gcta_scientific_summary
from postgwas.modules.gcta_gene.results import normalize_gcta_results
from postgwas.modules.gcta_gene.service import (
    _GCTA_MAXIMUM_SET_VARIANTS,
    _prepare_analysis_set_list,
    _resolved_configuration,
    run_gcta_gene_direct,
)
from postgwas.pipeline.runners import run_gcta_gene_runner


def _write_reference(tmp_path, variant="rs1"):
    prefix = tmp_path / "reference"
    Path(str(prefix) + ".bed").write_bytes(b"BED")
    Path(str(prefix) + ".fam").write_text("F1 I1 0 0 0 -9\n", encoding="utf-8")
    Path(str(prefix) + ".bim").write_text(
        "1\t%s\t0\t100\tG\tA\n" % variant, encoding="utf-8",
    )
    return prefix


def _write_gene_list(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("1\t50\t150\tENSG000001\n", encoding="utf-8")
    return path


def _write_set_list(tmp_path, variant="rs1"):
    path = tmp_path / "sets.txt"
    path.write_text("SET_1\n%s\nEND\n" % variant, encoding="utf-8")
    return path


def _write_input(tmp_path, method="fastbat_gene", variant="rs1"):
    path = tmp_path / "input.ma"
    path.write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "%s\tG\tA\t0.2\t0.1\t0.05\t0.01\t1000\n" % variant,
        encoding="utf-8",
    )
    return path


def _write_fake_gcta(tmp_path, version="1.94.1", version_exit_code=0):
    path = tmp_path / "fake_gcta"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        " print('GCTA v%s')\n"
        " raise SystemExit(%d)\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "if '--fastBAT' in sys.argv:\n"
        " if '--fastBAT-seg' in sys.argv:\n"
        "  Path(prefix + '.seg.fastbat').write_text("
        "'Chr Start End No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ '1 1 100000 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        " elif '--fastBAT-set-list' in sys.argv:\n"
        "  Path(prefix + '.fastbat').write_text("
        "'Set No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ 'SET_1 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        " else:\n"
        "  Path(prefix + '.gene.fastbat').write_text("
        "'Gene Chr Start End No.SNPs SNP_start SNP_end Chisq(Obs) Pvalue TopSNP.Pvalue TopSNP\\n'"
        "+ 'ENSG000001 1 50 150 1 rs1 rs1 4.0 0.01 0.01 rs1\\n')\n"
        "else:\n"
        " Path(prefix + '.gene.assoc.mbat').write_text("
        "'Gene Chr Start End No.SNPs SNP_start SNP_end TopSNP TopSNP_Pvalue No.Eigenvalues Chisq_mBAT P_mBATcombo P_mBAT Chisq_fastBAT P_fastBAT\\n'"
        "+ 'ENSG000001 1 50 150 1 rs1 rs1 rs1 0.01 1 4.0 0.02 0.03 4.0 0.01\\n')\n"
        % (version, version_exit_code),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_fake_bcftools(
    tmp_path, *, rows=(("1", 100, "rs1", "A", "G"),),
):
    path = tmp_path / "fake_bcftools"
    payload = "\n".join("\t".join(map(str, row)) for row in rows)
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '-l' in sys.argv:\n"
        " print('STUDY')\n"
        "else:\n"
        " print(%r)\n" % payload,
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_config(tmp_path, method, reference, genes, sets):
    path = tmp_path / "gcta_gene.yaml"
    path.write_text(
        "enabled: true\n"
        "method: %s\n"
        "genome_build: GRCh37\n"
        "reference:\n"
        "  prefix: '%s'\n"
        "  population: EUR\n"
        "gene_annotation:\n"
        "  file: '%s'\n"
        "set_annotation:\n"
        "  file: '%s'\n"
        % (method, reference, genes, sets),
        encoding="utf-8",
    )
    return path


def _args(tmp_path, method="fastbat_gene", *, version="1.94.1", variant="rs1"):
    reference = _write_reference(tmp_path)
    genes = _write_gene_list(tmp_path)
    sets = _write_set_list(tmp_path, variant)
    return Namespace(
        gcta_input_file=str(_write_input(tmp_path, method, variant)),
        dataset_id="STUDY",
        output_directory=str(tmp_path / "out"),
        run_config=str(_write_config(tmp_path, method, reference, genes, sets)),
        gcta=str(_write_fake_gcta(tmp_path, version)),
        threads=2,
        dry_run=False,
        resume=False,
        overwrite=False,
    )


def test_help_exposes_scientific_compatibility_and_configuration():
    parser = build_parser()
    help_text = parser.format_help()
    normalized_help = " ".join(help_text.split())
    assert (
        "--method {fastbat_gene,fastbat_segment,fastbat_set,mbat_combo}"
        in help_text
    )
    assert "--fastbat-set-list" in help_text
    assert "--gmt" in help_text
    assert "mutually exclusive with --fastbat-set-list" in normalized_help
    assert "--fastbat-segment-size-kb" in help_text
    build_actions = [
        action for action in parser._actions if "--genome-build" in action.option_strings
    ]
    assert len(build_actions) == 1
    assert "all coordinate-bearing inputs must match" in normalized_help
    assert "--gcta-reference-genome-build" not in help_text
    assert "--gene-list-genome-build" not in help_text
    assert "--gcta-reference-population" in help_text
    assert "postgwas config export" in help_text


def test_configurable_cli_values_do_not_own_defaults():
    parser = build_parser()
    for destination in (
        "gcta_input_file", "gcta_gene_method", "gcta_reference_prefix",
        "gene_list", "fastbat_set_list", "gmt", "genome_build",
        "gcta_reference_population", "gene_window_kb", "fastbat_segment_size_kb",
        "gmt_chromosome_label_policy", "gmt_duplicate_gene_policy",
        "gmt_unmapped_gene_policy", "gmt_empty_pathway_policy",
        "fastbat_oversized_set_policy",
        "gcta_reference_maf_min", "fastbat_ld_cutoff", "mbat_svd_gamma",
        "frequency_difference_max", "print_component_p_values", "write_snpset",
        "gcta_top_results", "gcta_nominal_alpha", "gcta_reporting_alpha",
        "gcta_fdr_alpha", "gcta_p_value_digits",
        "gcta", "bcftools", "vcf", "run_config", "resume", "overwrite",
        "dry_run", "dataset_id",
        "output_directory", "threads", "memory_gb", "seed",
        "gcta_coordinate_fallback", "gcta_chromosome_label_policy",
        "gcta_minimum_reference_overlap", "gcta_allow_strand_complement",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS


def test_cli_help_defaults_are_read_from_canonical_configuration(monkeypatch):
    defaults = load_configuration()
    module = defaults.modules.gcta_gene
    module.method = "mbat_combo"
    module.gene_window_kb = 73
    module.fastbat_ld_cutoff = 0.83
    module.print_component_p_values = False
    module.reporting.top_result_count = 7
    module.reporting.nominal_alpha = 0.02
    module.reporting.familywise_alpha = 0.01
    module.reporting.fdr_alpha = 0.03
    module.reporting.p_value_significant_digits = 6
    module.set_annotation.conversion.chromosome_label_policy = "strip_chr_prefix"
    defaults.resources.executables.gcta = "configured-gcta"
    monkeypatch.setattr(gcta_gene_cli, "load_configuration", lambda: defaults)

    parser = gcta_gene_cli.get_gcta_gene_parser(direct_controls=True)
    expected_by_destination = {
        "gcta_gene_method": "mbat_combo",
        "gene_window_kb": 73,
        "fastbat_ld_cutoff": 0.83,
        "print_component_p_values": False,
        "gcta_top_results": 7,
        "gcta_nominal_alpha": 0.02,
        "gcta_reporting_alpha": 0.01,
        "gcta_fdr_alpha": 0.03,
        "gcta_p_value_digits": 6,
        "gmt_chromosome_label_policy": "strip_chr_prefix",
        "fastbat_oversized_set_policy": "omit",
        "gcta": "configured-gcta",
    }
    for destination, expected in expected_by_destination.items():
        action = next(item for item in parser._actions if item.dest == destination)
        assert format_cli_default(
            expected, label="Configured YAML default",
        ) in action.help


def test_output_help_uses_canonical_run_defaults_without_argparse_defaults():
    defaults = load_configuration()
    parser = build_parser()
    expected_by_destination = {
        "dataset_id": defaults.run.dataset_id,
        "output_directory": defaults.run.output_directory,
    }
    for destination, expected in expected_by_destination.items():
        action = next(item for item in parser._actions if item.dest == destination)
        assert format_cli_default(
            expected, label="Configured YAML default",
        ) in action.help
        assert action.default == argparse.SUPPRESS


def test_pipeline_runner_defers_complete_cli_validation_to_service(
    tmp_path, monkeypatch,
):
    context = {"formatter": {"gcta_gene": {}}}
    args = Namespace(
        output_directory=str(tmp_path / "pipeline"),
        _step_num="02",
        gcta_gene_method="fastbat_set",
    )
    captured = {}
    expected = object()

    def reject_partial_configuration(*_args, **_kwargs):
        raise AssertionError("runner performed partial GCTA configuration validation")

    def fake_service(service_args, service_context):
        captured["has_input_override"] = hasattr(service_args, "gcta_input_file")
        captured["output_directory"] = service_args.output_directory
        captured["context"] = service_context
        return expected

    monkeypatch.setattr(
        "postgwas.pipeline.runners.load_module_configuration",
        reject_partial_configuration,
    )
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_gene_direct",
        fake_service,
    )

    result = run_gcta_gene_runner(args, context)

    assert result is expected
    assert captured == {
        "has_input_override": False,
        "output_directory": str(tmp_path / "pipeline" / "02_gcta_gene"),
        "context": context,
    }
    assert args.output_directory == str(tmp_path / "pipeline")


def test_pipeline_cli_path_values_resolve_as_canonical_configuration(tmp_path):
    reference = _write_reference(tmp_path)
    gene_list = _write_gene_list(tmp_path)
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    output = tmp_path / "pipeline"
    args = Namespace(
        gcta_gene_method="fastbat_set",
        gcta_reference_prefix=str(reference),
        gene_list=gene_list,
        gmt=gmt,
        genome_build="GRCh37",
        gcta_reference_population="EUR",
        dataset_id="STUDY",
        output_directory=output,
    )

    configuration = _resolved_configuration(args)
    module = configuration.modules.gcta_gene

    assert module.method == "fastbat_set"
    assert module.gene_annotation.file == str(gene_list)
    assert module.set_annotation.gmt_file == str(gmt)
    assert configuration.run.dataset_id == "STUDY"
    assert configuration.run.output_directory == output


@pytest.mark.parametrize(
    "method",
    [
        "fastbat_set",
        "mbat_combo",
    ],
)
def test_pipeline_service_selects_shared_formatted_input(tmp_path, method):
    args = _args(tmp_path, method)
    formatted_input = args.gcta_input_file
    del args.gcta_input_file
    context = RunContext({
        "formatter": {
            "gcta_gene": {"summary_statistics_input_file": formatted_input},
        },
    })

    result = run_gcta_gene_direct(args, context)

    assert result.metrics["method"] == method
    assert context["gcta_gene"] is result


def test_fastbat_command_uses_only_official_gene_test_options():
    module = load_module_configuration("gcta_gene")
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", "genes.txt", "out/STUDY", 4, module,
    )
    assert command == [
        "gcta64", "--bfile", "reference", "--maf", "0.01", "--thread-num", "4",
        "--out", "out/STUDY", "--fastBAT", "study.tsv", "--fastBAT-ld-cutoff",
        "0.9", "--fastBAT-gene-list", "genes.txt", "--fastBAT-wind", "50",
    ]
    assert "--mBAT-combo" not in command


def test_fastbat_segment_command_uses_official_segment_option(tmp_path):
    config = tmp_path / "segment.yaml"
    config.write_text("method: fastbat_segment\n", encoding="utf-8")
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", None, "out/STUDY", 4, module,
    )
    assert "--fastBAT-seg" in command
    assert command[command.index("--fastBAT-seg") + 1] == "100"
    assert "--fastBAT-gene-list" not in command
    assert "--fastBAT-set-list" not in command


def test_fastbat_set_command_uses_official_set_and_snpset_options(tmp_path):
    config = tmp_path / "set.yaml"
    config.write_text(
        "method: fastbat_set\nwrite_snpset: true\n", encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.tsv", "reference", "sets.txt", "out/STUDY", 4,
        module,
    )
    assert "--fastBAT-set-list" in command
    assert "--fastBAT-write-snpset" in command
    assert "--fastBAT-gene-list" not in command


def test_mbat_combo_command_includes_every_documented_method_option(tmp_path):
    config = tmp_path / "mbat.yaml"
    config.write_text(
        "method: mbat_combo\nwrite_snpset: true\n", encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene", config)
    command = build_gcta_command(
        "gcta64", "study.ma", "reference", "genes.txt", "out/STUDY", 4,
        module,
    )
    assert command == [
        "gcta64", "--bfile", "reference", "--maf", "0.01", "--thread-num", "4",
        "--out", "out/STUDY", "--mBAT-combo", "study.ma", "--mBAT-gene-list",
        "genes.txt", "--mBAT-wind", "50", "--mBAT-svd-gamma", "0.9",
        "--diff-freq", "0.2", "--fastBAT-ld-cutoff", "0.9",
        "--mBAT-print-all-p", "--mBAT-write-snpset",
    ]


def test_gcta_input_accepts_distinct_nonempty_indel_alleles(tmp_path):
    args = _args(tmp_path, "fastbat_gene")
    Path(args.gcta_input_file).write_text(
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\n"
        "rs1\tCAAATAAAT\tC\t0.2\t0.1\t0.05\t0.01\t1000\n",
        encoding="utf-8",
    )
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\trs1\t0\t100\tC\tCAAATAAAT\n", encoding="utf-8",
    )

    result = run_gcta_gene_direct(args)

    assert result.metrics["compatible_allele_pairs"] == 1


def test_mbat_combo_runs_once_and_reports_internal_fastbat(tmp_path, capsys):
    args = _args(tmp_path, "mbat_combo")
    result = run_gcta_gene_direct(args)

    assert result.metrics["method"] == "mbat_combo"
    assert result.metrics["tested_units"] == 1
    assert result.metrics["minimum_p_value"] == pytest.approx(0.02)
    top = result.metrics["scientific_summary"]["top_associations"][0]
    assert top["p_value"] == pytest.approx(0.02)
    assert top["component_p_values"] == {
        "P_mBAT": pytest.approx(0.03),
        "P_fastBAT": pytest.approx(0.01),
    }
    raw = result.artifacts["raw_results"].path
    normalized = result.artifacts["normalized_results"].path
    assert raw.name == "STUDY.gene.assoc.mbat" and raw.is_file()
    assert "P_fastBAT" in normalized.read_text(encoding="utf-8").splitlines()[0]
    log = tmp_path / "out" / "logs" / "STUDY_mbat_combo_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    assert log_text.count("--mBAT-combo") == 1
    assert "--mBAT-print-all-p" in log_text
    assert "--fastBAT," not in log_text
    assert "gcta_scientific_summary" in log_text
    assert "gcta_top_association" in log_text
    terminal = capsys.readouterr().out
    assert "GCTA mBAT-combo gene-association summary" in terminal
    assert "P_mBATcombo 0.02" in terminal
    assert "P_fdr_bh 0.02" in terminal
    assert "P_bonferroni 0.02" in terminal
    assert "P_mBAT 0.03" in terminal
    assert "P_fastBAT 0.01" in terminal
    assert not list((tmp_path / "out" / ".partial").glob("*"))
    resolved = yaml.safe_load(
        (tmp_path / "out" / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert list(resolved["modules"]) == ["formatting", "gcta_gene"]
    assert set(resolved["resources"]) == {
        "executables", "genomes", "populations",
    }
    assert list(resolved["resources"]["genomes"]) == ["GRCh37"]
    assert list(resolved["resources"]["populations"]) == ["EUR"]


def test_fastbat_run_validates_and_normalizes_results(tmp_path):
    result = run_gcta_gene_direct(_args(tmp_path))
    assert result.metrics["gcta_version"] == "1.94.1"
    assert result.metrics["overlapping_variants"] == 1
    assert result.metrics["tested_units"] == 1
    assert result.artifacts["raw_results"].path.name == "STUDY.gene.fastbat"
    assert result.artifacts["normalized_results"].path.is_file()
    assert result.artifacts["completion_manifest"].path.is_file()
    summary = result.metrics["scientific_summary"]
    assert summary["tested_units"] == 1
    assert summary["bonferroni_threshold"] == pytest.approx(0.05)
    assert summary["nominal_significant"] == 1
    assert summary["fdr_bh_significant"] == 1
    assert summary["bonferroni_significant"] == 1
    corrected = pl.read_csv(
        result.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert corrected["P_bonferroni"].to_list() == pytest.approx([0.01])
    assert corrected["P_fdr_bh"].to_list() == pytest.approx([0.01])
    assert corrected["significant_nominal"].to_list() == [True]
    assert corrected["significant_bonferroni"].to_list() == [True]
    assert corrected["significant_fdr_bh"].to_list() == [True]


def test_scientific_summary_ranks_results_and_counts_bonferroni_hits(tmp_path):
    raw_result = tmp_path / "raw.tsv"
    raw_result.write_text(
        "Gene\tChr\tStart\tEnd\tNo.SNPs\tPvalue\n"
        "GENE_C\t1\t300\t350\t3\t0.2\n"
        "GENE_A\t1\t100\t150\t2\t0.001\n"
        "GENE_B\t1\t200\t250\t4\t0.01\n",
        encoding="utf-8",
    )
    module = load_module_configuration("gcta_gene")
    module.reporting.top_result_count = 2
    result = tmp_path / "results.tsv"
    normalize_gcta_results(raw_result, result, module)

    summary = build_gcta_scientific_summary(
        result,
        module,
        {
            "tested_units": 3,
            "variants": 100,
            "overlapping_variants": 90,
            "unresolved_variants_removed": 10,
            "genes": 3,
            "shared_chromosomes": ["1", "2"],
        },
    )

    assert [item["identifier"] for item in summary["top_associations"]] == [
        "GENE_A", "GENE_B",
    ]
    assert summary["bonferroni_threshold"] == pytest.approx(0.05 / 3)
    assert summary["nominal_significant"] == 2
    assert summary["fdr_bh_significant"] == 2
    assert summary["bonferroni_significant"] == 2
    assert summary["overlap_fraction"] == pytest.approx(0.9)
    assert summary["source_units"] == 3
    assert summary["result_chromosomes"] == ["1"]
    assert summary["missing_result_chromosomes"] == ["2"]
    assert "10 GWAS variants" in summary["warnings"][0]
    assert "shared chromosome(s): 2" in summary["warnings"][1]
    corrected = pl.read_csv(result, separator="\t")
    assert corrected["P_bonferroni"].to_list() == pytest.approx([
        0.6, 0.003, 0.03,
    ])
    assert corrected["P_fdr_bh"].to_list() == pytest.approx([
        0.2, 0.003, 0.015,
    ])


def test_completed_gcta_result_resumes_only_after_manifest_validation(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path)
    first = run_gcta_gene_direct(args)
    capsys.readouterr()
    args.resume = True
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("validated result was rerun"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["tested_units"] == first.metrics["tested_units"]
    assert resumed.metrics["resumed"] is True
    assert (
        "GCTA outputs validated; continuing from completed step"
        in capsys.readouterr().out
    )
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "reason=validated_resume" in log.read_text(encoding="utf-8")


def test_reporting_changes_do_not_invalidate_completed_analysis(
    tmp_path, monkeypatch, capsys,
):
    args = _args(tmp_path)
    run_gcta_gene_direct(args)
    capsys.readouterr()
    args.resume = True
    args.gcta_top_results = 1
    args.gcta_nominal_alpha = 0.005
    args.gcta_reporting_alpha = 0.005
    args.gcta_fdr_alpha = 0.005
    args.gcta_p_value_digits = 6
    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.run_gcta_command",
        lambda *_args, **_kwargs: pytest.fail("reporting change reran GCTA"),
    )

    resumed = run_gcta_gene_direct(args)

    assert resumed.metrics["resumed"] is True
    summary = resumed.metrics["scientific_summary"]
    assert summary["nominal_alpha"] == 0.005
    assert summary["familywise_alpha"] == 0.005
    assert summary["fdr_alpha"] == 0.005
    assert summary["nominal_significant"] == 0
    assert summary["bonferroni_significant"] == 0
    assert summary["fdr_bh_significant"] == 0
    corrected = pl.read_csv(
        resumed.artifacts["normalized_results"].path,
        separator="\t",
    )
    assert corrected["significant_nominal"].to_list() == [False]
    assert corrected["significant_bonferroni"].to_list() == [False]
    assert corrected["significant_fdr_bh"].to_list() == [False]
    terminal = capsys.readouterr().out
    assert "GCTA outputs validated; continuing from completed step" in terminal
    assert "GCTA fastBAT gene-association summary" in terminal
    assert "Nominally significant associations" in terminal
    assert "FDR-significant associations" in terminal
    assert "Bonferroni threshold" in terminal
    assert "Bonferroni-significant associations" in terminal
    assert "Raw p ≤" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "gcta_multiple_testing_corrections" in log.read_text(encoding="utf-8")


def test_completed_gcta_result_rejects_changed_output_on_resume(tmp_path):
    args = _args(tmp_path)
    result = run_gcta_gene_direct(args)
    result.artifacts["raw_results"].path.write_text(
        "changed\n", encoding="utf-8",
    )
    args.resume = True

    with pytest.raises(GctaGeneError, match="completion outputs changed"):
        run_gcta_gene_direct(args)


def test_pipeline_vcf_reconciles_rsid_input_to_exact_bim_identifier(tmp_path):
    args = _args(tmp_path, "fastbat_set", variant="rs1")
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\t1_100_A_G\t0\t100\tG\tA\n", encoding="utf-8",
    )
    (tmp_path / "sets.txt").write_text(
        "SET_1\n1_100_A_G\nEND\n", encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    args.bcftools = str(_write_fake_bcftools(tmp_path))

    result = run_gcta_gene_direct(args)

    assert result.metrics["direct_id_matches"] == 0
    assert result.metrics["coordinate_and_allele_matches"] == 1
    assert result.metrics["variant_ids_replaced"] == 1
    harmonised = result.artifacts["harmonised_input"].path
    assert harmonised.read_text(encoding="utf-8").splitlines() == [
        "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN",
        "1_100_A_G\tG\tA\t0.2\t0.1\t0.05\t0.01\t1000",
    ]
    assert result.artifacts["analysis_set_list"].path.read_text(
        encoding="utf-8"
    ).split() == ["SET_1", "1_100_A_G", "END"]


def test_coordinate_reconciliation_enforces_configured_minimum_overlap(tmp_path):
    args = _args(tmp_path)
    input_path = Path(args.gcta_input_file)
    input_path.write_text(
        input_path.read_text(encoding="utf-8")
        + "rs2\tC\tT\t0.3\t0.2\t0.1\t0.02\t1000\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    args.bcftools = str(_write_fake_bcftools(
        tmp_path,
        rows=(("1", 100, "rs1", "A", "G"), ("1", 200, "rs2", "C", "T")),
    ))
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["variant_harmonisation"] = {"minimum_overlap_fraction": 0.75}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(GctaGeneError, match="Only 1/2 unique GCTA input variants"):
        run_gcta_gene_direct(args)


def test_direct_id_match_rejects_vcf_reference_coordinate_conflict(tmp_path):
    args = _args(tmp_path)
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"VCF")
    args.vcf = str(vcf)
    args.bcftools = str(_write_fake_bcftools(
        tmp_path, rows=(("1", 101, "rs1", "A", "G"),),
    ))

    with pytest.raises(GctaGeneError, match="direct ID matches have different"):
        run_gcta_gene_direct(args)


def test_old_gcta_version_fails_and_is_logged(tmp_path):
    args = _args(tmp_path, version="1.93.3")
    with pytest.raises(GctaGeneError, match="too old"):
        run_gcta_gene_direct(args)
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "FAILED" in log.read_text(encoding="utf-8")


def test_valid_version_banner_is_accepted_when_probe_exits_nonzero(tmp_path):
    args = _args(tmp_path)
    args.gcta = str(_write_fake_gcta(tmp_path, version_exit_code=1))

    result = run_gcta_gene_direct(args)

    assert result.metrics["gcta_version"] == "1.94.1"
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "probe_exit_code=1" in log.read_text(encoding="utf-8")


def test_zero_reference_overlap_fails_before_analysis(tmp_path):
    args = _args(tmp_path, variant="rs_missing")
    with pytest.raises(GctaGeneError, match="No GCTA input variant identifiers"):
        run_gcta_gene_direct(args)
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_gene_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")


def test_mbat_combo_rejects_any_incompatible_reference_allele_pair(tmp_path):
    args = _args(tmp_path, "mbat_combo")
    input_path = Path(args.gcta_input_file)
    input_path.write_text(
        input_path.read_text(encoding="utf-8")
        + "rs2\tC\tA\t0.3\t0.2\t0.1\t0.02\t1000\n",
        encoding="utf-8",
    )
    reference = Path(str(tmp_path / "reference") + ".bim")
    reference.write_text(
        reference.read_text(encoding="utf-8")
        + "1\trs2\t0\t110\tG\tA\n",
        encoding="utf-8",
    )

    with pytest.raises(GctaGeneError, match="1 overlapping.*cannot be matched"):
        run_gcta_gene_direct(args)


def test_duplicate_nested_genome_build_declarations_are_rejected(tmp_path):
    config = tmp_path / "invalid.yaml"
    config.write_text(
        "genome_build: GRCh37\nreference:\n  genome_build: GRCh37\n",
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="Extra inputs are not permitted"):
        load_module_configuration("gcta_gene", config)


@pytest.mark.parametrize(
    "text",
    (
        "fastbat_ld_cutoff: 1.1\n",
        "mbat_svd_gamma: 0\n",
        "frequency_difference_max: .nan\n",
        "variant_harmonisation:\n  minimum_overlap_fraction: 0\n",
        "reporting:\n  top_result_count: 0\n",
        "reporting:\n  nominal_alpha: 0\n",
        "reporting:\n  familywise_alpha: 1.1\n",
        "reporting:\n  fdr_alpha: .nan\n",
        "reporting:\n  p_value_significant_digits: 11\n",
        "reporting:\n  correction_columns:\n    fdr_bh_adjusted_p: P_bonferroni\n",
        "reporting:\n  chromosome_columns:\n    fastbat_gene: Missing\n",
    ),
)
def test_scientific_ranges_are_schema_validated(tmp_path, text):
    config = tmp_path / "invalid.yaml"
    config.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigurationError):
        load_module_configuration("gcta_gene", config)


def test_reference_and_gene_column_orders_are_configuration_driven(tmp_path):
    args = _args(tmp_path, "mbat_combo")
    reference = Path(str(tmp_path / "reference") + ".bim")
    reference.write_text("rs1\t1\tA\tG\t100\t0\n", encoding="utf-8")
    genes = tmp_path / "genes.txt"
    genes.write_text("ENSG000001\t150\t1\t50\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["reference"]["bim_columns"] = [
        "variant_id", "chromosome", "allele2", "allele1", "position",
        "genetic_distance",
    ]
    config["gene_annotation"]["columns"] = ["gene", "end", "chromosome", "start"]
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    result = run_gcta_gene_direct(args)

    assert result.metrics["compatible_allele_pairs"] == 1
    assert result.metrics["genes"] == 1


@pytest.mark.parametrize(
    ("method", "result_name", "unit_label"),
    (
        ("fastbat_segment", "STUDY.seg.fastbat", "segments"),
        ("fastbat_set", "STUDY.fastbat", "sets"),
    ),
)
def test_all_fastbat_analysis_modes_run_and_normalize(
    tmp_path, method, result_name, unit_label, capsys,
):
    result = run_gcta_gene_direct(_args(tmp_path, method))

    assert result.artifacts["raw_results"].path.name == result_name
    assert result.metrics["tested_units"] == 1
    assert result.metrics["unit_label"] == unit_label
    assert result.artifacts["normalized_results"].path.is_file()
    assert "Top associated %s" % unit_label in capsys.readouterr().out


def test_fastbat_set_rejects_an_unterminated_set_before_execution(tmp_path):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text("SET_1\nrs1\n", encoding="utf-8")

    with pytest.raises(GctaGeneError, match="missing its terminating END"):
        run_gcta_gene_direct(args)

    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")


def test_fastbat_set_converts_gmt_with_exact_bim_identifiers(tmp_path):
    args = _args(tmp_path, "fastbat_set", variant="custom_1_100_A_G")
    Path(str(tmp_path / "reference") + ".bim").write_text(
        "1\tcustom_1_100_A_G\t0\t100\tG\tA\n", encoding="utf-8",
    )
    genes = tmp_path / "genes.txt"
    genes.write_text("1\t50\t150\tGENE1\n", encoding="utf-8")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tGENE1\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    result = run_gcta_gene_direct(args)

    prepared = tmp_path / "out" / "prepared_sets" / "STUDY"
    set_text = (prepared / "pathways.fastbat.set").read_text(encoding="utf-8")
    assert set_text.split() == ["SET_1", "custom_1_100_A_G", "END"]
    assert result.artifacts["prepared_set_resource"].path == prepared
    assert result.artifacts["analysis_set_list"].path == (
        prepared / "analysis.fastbat.set"
    )
    assert result.metrics["matched_set_variants"] == 1
    manifest = yaml.safe_load(
        (prepared / "resource_manifest.yaml").read_text(encoding="utf-8")
    )
    assert manifest["resource"]["genome_build"] == "GRCh37"
    assert manifest["policies"]["variant_identifier_policy"].startswith("copy PLINK BIM")


def test_fastbat_set_omits_zero_overlap_sets_from_analysis_input(tmp_path, capsys):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text(
        "EMPTY_SET\nrs_missing\nEND\n\nSET_1\nrs1\nEND\n",
        encoding="utf-8",
    )

    result = run_gcta_gene_direct(args)

    analysis_set = result.artifacts["analysis_set_list"].path
    assert analysis_set.read_text(encoding="utf-8").split() == [
        "SET_1", "rs1", "END",
    ]
    assert result.metrics["input_sets"] == 2
    assert result.metrics["sets"] == 1
    assert result.metrics["omitted_empty_sets"] == 1
    assert result.metrics["omitted_empty_set_examples"] == ["EMPTY_SET"]
    assert len(result.metrics["scientific_summary"]["warnings"]) == 1
    terminal = capsys.readouterr().out
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in terminal
    assert "1 custom sets had no variants shared" in terminal
    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    log_text = log.read_text(encoding="utf-8")
    assert "fastbat_set_intersection" in log_text
    assert "1 custom sets had no variants shared" in log_text


def test_fastbat_set_empty_overlap_error_policy_stops_before_gcta(tmp_path):
    args = _args(tmp_path, "fastbat_set")
    (tmp_path / "sets.txt").write_text(
        "EMPTY_SET\nrs_missing\nEND\n\nSET_1\nrs1\nEND\n",
        encoding="utf-8",
    )
    config_path = Path(args.run_config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["set_annotation"]["conversion"] = {"empty_pathway_policy": "error"}
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(GctaGeneError, match="1 fastBAT sets have no variants"):
        run_gcta_gene_direct(args)

    log = tmp_path / "out" / "logs" / "STUDY_fastbat_set_gcta_gene.log"
    assert "--fastBAT" not in log.read_text(encoding="utf-8")


@pytest.mark.parametrize("policy", ["omit", "error"])
def test_fastbat_set_enforces_gcta_variant_limit_before_execution(tmp_path, policy):
    variants = ["rs%d" % index for index in range(_GCTA_MAXIMUM_SET_VARIANTS + 1)]
    source = tmp_path / "sets.txt"
    source.write_text(
        "OVERSIZED\n%s\nEND\n\nVALID\nrs_valid\nEND\n"
        % "\n".join(variants),
        encoding="utf-8",
    )
    identifiers = set(variants) | {"rs_valid"}
    destination = tmp_path / "analysis.set"

    if policy == "error":
        with pytest.raises(GctaGeneError, match="20,000-variant hard limit"):
            _prepare_analysis_set_list(
                source,
                destination,
                identifiers,
                identifiers,
                "omit",
                policy,
            )
        return

    _, metrics = _prepare_analysis_set_list(
        source,
        destination,
        identifiers,
        identifiers,
        "omit",
        policy,
    )
    assert destination.read_text(encoding="utf-8").split() == [
        "VALID", "rs_valid", "END",
    ]
    assert metrics["sets"] == 1
    assert metrics["omitted_oversized_sets"] == 1
    assert metrics["omitted_oversized_set_examples"] == [
        "OVERSIZED (20001)"
    ]
    assert metrics["analysis_set_variants"] == 1


def test_fastbat_set_rejects_gmt_and_prepared_set_together(tmp_path):
    args = _args(tmp_path, "fastbat_set")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"]["gmt_file"] = str(gmt)
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ConfigurationError, match="mutually exclusive"):
        run_gcta_gene_direct(args)


def test_prepared_set_input_bypasses_gmt_conversion(tmp_path, monkeypatch):
    args = _args(tmp_path, "fastbat_set")
    gmt = tmp_path / "configured_pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")
    args.fastbat_set_list = str(tmp_path / "sets.txt")

    def unexpected_conversion(_):
        raise AssertionError("GMT converter must not run for --fastbat-set-list")

    monkeypatch.setattr(
        "postgwas.modules.gcta_gene.service.prepare_resource",
        unexpected_conversion,
    )
    result = run_gcta_gene_direct(args)

    assert result.metrics["tested_units"] == 1
    assert "prepared_set_resource" not in result.artifacts


def test_gmt_is_rejected_for_non_set_method(tmp_path):
    args = _args(tmp_path, "fastbat_gene")
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("SET_1\tdescription\tENSG000001\n", encoding="utf-8")
    config = yaml.safe_load(Path(args.run_config).read_text(encoding="utf-8"))
    config["set_annotation"] = {"file": None, "gmt_file": str(gmt)}
    Path(args.run_config).write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(GctaGeneError, match="valid only with --method fastbat_set"):
        run_gcta_gene_direct(args)


def test_dry_run_removes_empty_staging_directory(tmp_path):
    args = _args(tmp_path)
    args.dry_run = True
    result = run_gcta_gene_direct(args)
    assert result.metrics["dry_run"] is True
    partial = tmp_path / "out" / ".partial" / "STUDY_fastbat_gene"
    assert not partial.exists()
