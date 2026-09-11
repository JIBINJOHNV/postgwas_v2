import argparse
from argparse import Namespace
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from postgwas.cli.common import get_genome_build_parser
from postgwas.config import load_configuration, load_module_configuration
from postgwas.core.contracts import RunContext
from postgwas.core.errors import ConfigurationError
from postgwas.core.input_validation import InputValidationSession
from postgwas.core.validation_reporting import FileValidationDisplay
from postgwas.modules.mixer.cli import build_parser
from postgwas.modules.mixer.service import (
    MixerError,
    MixerPipelineResources,
    preflight_mixer_pipeline,
    run_mixer_direct,
    validate_reference_pattern,
)
from postgwas.pipeline.runners import run_mixer_runner
from preflight_support import pipeline_input_vcf_evidence


def _module_config(tmp_path, *, analysis="univariate", extra=""):
    path = tmp_path / "mixer.yaml"
    path.write_text(
        "enabled: true\n"
        "analysis: %s\n"
        "genome_build: GRCh37\n"
        "chromosomes: ['1']\n"
        "%s" % (analysis, extra),
        encoding="utf-8",
    )
    return path


def _mixer_input(tmp_path):
    path = tmp_path / "study.sumstats.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tCHR\tBP\tA1\tA2\tN\tZ\nrs1\t1\t1\tA\tG\t1000\t2\n")
    return path


def test_mixer_chromosome_resources_have_one_compact_validation_card(
    tmp_path, capsys,
):
    pattern = str(tmp_path / "reference_chr@.bim")
    paths = []
    for chromosome in ("1", "2", "3"):
        path = Path(pattern.replace("@", chromosome))
        path.write_text("reference\n", encoding="utf-8")
        paths.append(path.resolve())

    with InputValidationSession() as session:
        display = FileValidationDisplay(session, load_configuration())
        observed = validate_reference_pattern(
            pattern, "BIM", ("1", "2", "3"), "@",
        )
        display.flush()

    assert observed == tuple(paths)
    screen = capsys.readouterr().out
    assert "MiXeR chromosome resources — AVAILABLE — availability only" in screen
    assert "Resource type" in screen and "BIM" in screen
    assert "Chromosome files" in screen and "3 / 3" in screen
    assert all(path.name not in screen for path in paths)


def _references(tmp_path, *, loadlib=False):
    paths = {}
    for label, suffix in (
        ("bim", "bim"), ("ld", "ld"), ("annotation", "annot.gz"),
    ):
        path = tmp_path / ("reference_chr1.%s" % suffix)
        path.write_text("reference\n", encoding="utf-8")
        paths[label] = str(path).replace("chr1", "chr@")
    if loadlib:
        path = tmp_path / "reference_chr1.bin"
        path.write_text("load library\n", encoding="utf-8")
        paths["loadlib"] = str(path).replace("chr1", "chr@")
    return paths


def _go_files(tmp_path):
    paths = {}
    records = {
        "baseline": [("coding_genes", "GENE1")],
        "model": [("coding_genes", "GENE1"), ("GENE1", "GENE1")],
        "test": [("SET_A", "GENE1"), ("SET_B", "GENE1"), ("SET_C", "GENE1")],
    }
    for label, rows in records.items():
        path = tmp_path / ("%s.tsv" % label)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["GO", "GENE", "CHR", "FROM", "TO"])
            for identifier, gene in rows:
                writer.writerow([identifier, gene, 1, 1, 2])
        paths[label] = str(path)
    return paths


def _base_args(tmp_path, output, config, references, **values):
    arguments = {
        "mixer_input_file": str(_mixer_input(tmp_path)),
        "dataset_id": "STUDY",
        "output_directory": str(output),
        "run_config": str(config),
        "bim_file_pattern": references["bim"],
        "ld_file_pattern": references.get("ld"),
        "threads": 1,
        "seed": 31,
        "dry_run": False,
        "resume": False,
        "overwrite": False,
    }
    arguments.update(values)
    return Namespace(**arguments)


def _canonical_log(output):
    logs = list((output / "logs").glob("*/*_mixer.log"))
    assert len(logs) == 1
    return logs[0]


def _write_fake_mixer(path):
    path.write_text(
        "import csv, json, sys\n"
        "from pathlib import Path\n"
        "command = sys.argv[1]\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "if command == 'split_sumstats':\n"
        " chrom = sys.argv[sys.argv.index('--chr2use') + 1]\n"
        " first, last = (chrom.split('-') + [chrom])[:2]\n"
        " for number in range(int(first), int(last) + 1):\n"
        "  target = Path(prefix.replace('@', str(number)))\n"
        "  target.parent.mkdir(parents=True, exist_ok=True)\n"
        "  target.write_text('SNP\\tCHR\\tBP\\tA1\\tA2\\tN\\tZ\\nrs1\\t1\\t1\\tA\\tG\\t1000\\t2\\n')\n"
        " sys.exit(0)\n"
        "seed = int(sys.argv[sys.argv.index('--seed') + 1])\n"
        "if command == 'plsa':\n"
        " base = '--gsa-base' in sys.argv\n"
        " options = {'gsa_base': base, 'gsa_full': not base, 'seed': seed, "
        "'num_snp': 100, 'num_tag': 80, 'trait1_nval': 1000, 'threads': [1]}\n"
        " Path(prefix + '.json').write_text(json.dumps({'analysis': 'plsa', 'options': options}))\n"
        " Path(prefix + '.log').write_text('MiXeR v2.2.1 software\\n')\n"
        " if base:\n"
        "  Path(prefix + '.snps.csv').write_text('snp\\nrs1\\n')\n"
        "  Path(prefix + '.weights').write_text('weights\\n')\n"
        " else:\n"
        "  columns = ['h2','se_h2','h2_base','enrich','se_enrich','snps','GO','NGENES',"
        "'loglike_diff','loglike_df','loglike_aic']\n"
        "  rows = [[.2,.01,.1,1.5,.2,10,'SET_A',2,3,1,4],"
        "[.2,.01,.1,9,.2,10,'SET_B',2,0,1,-1],"
        "[.2,.01,.1,3,.2,10,'SET_C',2,2,1,2]]\n"
        "  with open(prefix + '.go_test_enrich.csv', 'w', newline='') as handle:\n"
        "   writer = csv.writer(handle, delimiter='\\t'); writer.writerow(columns); writer.writerows(rows)\n"
        " sys.exit(0)\n"
        "ci = {'pi': {'point_estimate': .002}, 'nc': {'point_estimate': 20000.0}, "
        "'nc@p9': {'point_estimate': 7000.0}, 'sig2_beta': {'point_estimate': .00003}, "
        "'sig2_zero': {'point_estimate': 1.05 if command == 'fit1' else 1.06}, "
        "'h2': {'point_estimate': .18}}\n"
        "result = {'analysis': 'univariate', 'ci': ci, "
        "'options': {'seed': seed, 'num_snp': 10000000, 'num_tag': 1, 'threads': [1]}, "
        "'optimize': [['neldermead', {'success': True, 'AIC': 20.0, 'BIC': 30.0}]], "
        "'inft_optimize': [['infinitesimal', {'success': True, 'AIC': 40.0, 'BIC': 45.0}]]}\n"
        "if command == 'test1':\n"
        " result['qqplot'] = {'n_snps': 1, 'data_logpvec': [float('inf')], 'model_logpvec': [1.0]}\n"
        "Path(prefix + '.json').write_text(json.dumps(result))\n"
        "ld_pattern = sys.argv[sys.argv.index('--ld-file') + 1]\n"
        "ld_file = ld_pattern.replace('@', '1')\n"
        "with open(prefix + '.log', 'w') as handle:\n"
        " handle.write('MiXeR v2.2.1 software\\n')\n"
        " handle.write('<init(...); elapsed time 1ms\\n')\n"
        " handle.write('>load_ld_matrix(filename=%s)\\n' % ld_file)\n"
        " handle.write(\"--fit-sequence: ['diffevo-fast', 'neldermead']\\n\")\n"
        " handle.write('fit_type==diffevo-fast...\\n')\n"
        " handle.write('<calc_unified_univariate_cost_gaussian(...), cost=1, elapsed time 1ms\\n')\n"
        " handle.write('fit_type==diffevo-fast done (success=True)\\n')\n"
        " handle.write('Done\\n')\n"
        " if command == 'fit1':\n"
        "  handle.write('1 lines matched via CHR:BP:A1:A2 code (not SNP rs#)\\n')\n"
        "  handle.write('0 lines were ignored as RS# did not match reference file.\\n')\n"
        "  handle.write('0 variants were are strand-ambiguous, removing them from analysis\\n')\n"
        "  handle.write('0 variants had flipped A1/A2 alleles; sign of z-score was flipped.\\n')\n"
        "  handle.write('Found 1 variants with well-defined Z and N\\n')\n",
        encoding="utf-8",
    )


def _write_failing_mixer(path):
    path.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "ld_pattern = sys.argv[sys.argv.index('--ld-file') + 1]\n"
        "log = Path(prefix + '.log')\n"
        "log.parent.mkdir(parents=True, exist_ok=True)\n"
        "log.write_text(\n"
        " '<init(...); elapsed time 1ms\\n'\n"
        " + '>load_ld_matrix(filename=%s)\\n' % ld_pattern.replace('@', '1')\n"
        " + \"--fit-sequence: ['diffevo-fast', 'neldermead']\\n\"\n"
        " + 'fit_type==diffevo-fast...\\n'\n"
        " + '<calc_unified_univariate_cost_gaussian(...), cost=1, elapsed time 1ms\\n',\n"
        " encoding='utf-8',\n"
        ")\n"
        "sys.exit(2)\n",
        encoding="utf-8",
    )


def _write_fake_figures(path):
    path.write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "prefix = sys.argv[sys.argv.index('--out') + 1]\n"
        "for suffix in ('.csv', '.power.png', '.qq.png', '.qqbin.png'):\n"
        " Path(prefix + suffix).write_text(' '.join(sys.argv))\n",
        encoding="utf-8",
    )


def test_mixer_help_is_single_trait_and_exposes_gsa_resources():
    help_text = build_parser().format_help()
    assert help_text.startswith(
        "Usage: postgwas mixer --mixer-input-file PATH\n"
        "                      [--analysis {univariate,gsa,all}] [options]"
    )
    assert "--analysis {univariate,gsa,all}" in help_text
    assert "--gsa-annotation-file-pattern PATTERN" in help_text
    assert "--gsa-baseline-go-file PATH" in help_text
    assert "--gsa-model-go-file PATH" in help_text
    assert "--gsa-test-go-file PATH" in help_text
    assert "--threads N" in help_text
    assert "--dataset-id NAME" in help_text
    assert "Export a complete reloadable MiXeR configuration:" in help_text
    assert "  postgwas config export \\\n    --module mixer \\\n" in help_text
    assert "Set analysis: all in mixer.yaml" in help_text
    assert "GSA-MiXeR reference inputs:" in help_text
    description_markers = (
        "Total CPU-thread budget",
        "Per-chromosome SNP annotation",
        "GSA-MiXeR baseline annotation",
        "Single-trait analysis to run",
    )
    description_columns = {
        next(line.index(marker) for line in help_text.splitlines() if marker in line)
        for marker in description_markers
    }
    assert len(description_columns) == 1
    for forbidden in ("--trait2", "fit2", "test2", "bivariate", "--replicates", "--extract"):
        assert forbidden not in help_text.lower()


def test_pipeline_mixer_help_hides_only_formatter_created_input():
    completed = subprocess.run(
        [
            sys.executable, "-m", "postgwas", "pipeline", "--modules", "mixer",
            "--help",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    help_text = completed.stdout
    assert "Usage: postgwas pipeline [options] --modules mixer" in help_text
    assert "1) formatter" in help_text and "2) mixer" in help_text
    assert "Configuration:" in help_text
    assert "Workflow examples" in help_text
    assert "  postgwas config export \\\n    --pipeline mixer \\\n" in help_text
    assert "    --run-config analysis_pipeline.yaml" in help_text
    assert "--mixer-input-file" not in help_text
    for external_resource in (
        "--bim-file-pattern",
        "--ld-file-pattern",
        "GSA-MiXeR reference inputs:",
        "--analysis",
    ):
        assert external_resource in help_text
    description_markers = (
        "Total CPU-thread budget",
        "Required: Harmonised GWAS-VCF file",
        "Required: Short, unique name",
        "YAML file containing settings",
    )
    description_columns = {
        next(line.index(marker) for line in help_text.splitlines() if marker in line)
        for marker in description_markers
    }
    assert len(description_columns) == 1


def test_mixer_configurable_cli_actions_do_not_own_defaults():
    parser = build_parser()
    for destination in (
        "analysis", "genome_build", "bim_file_pattern", "ld_file_pattern",
        "gsa_annotation_file_pattern", "gsa_loadlib_file_pattern",
        "gsa_baseline_go_file", "gsa_model_go_file", "gsa_test_go_file",
        "mixer_backend", "mixer", "mixer_figures", "mixer_container_runtime",
        "mixer_container_image", "mixer_container_platform", "run_config",
        "resume", "overwrite", "dry_run", "dataset_id", "output_directory",
        "threads", "memory_gb", "seed",
    ):
        action = next(item for item in parser._actions if item.dest == destination)
        assert action.default == argparse.SUPPRESS


def test_omitted_genome_build_is_resolved_from_configuration():
    args = build_parser().parse_args([
        "--mixer-input-file", __file__, "--dataset-id", "STUDY",
        "--output-directory", ".",
    ])
    assert not hasattr(args, "genome_build")


def test_shared_genome_build_parser_uses_configured_names():
    parser = get_genome_build_parser(
        available_builds=["ReferenceA", "ReferenceB"],
        default_build="ReferenceB",
        suppress_default=True,
    )
    assert parser.parse_args(["--genome-build", "ReferenceA"]).genome_build == "ReferenceA"
    assert "Default: ReferenceB" in parser.format_help()


def test_pipeline_preflight_validates_mixer_resources_and_runner_reuses_them(
    tmp_path, monkeypatch,
):
    references = _references(tmp_path)
    output = tmp_path / "out"
    fake_mixer = tmp_path / "fake_mixer.py"
    _write_fake_mixer(fake_mixer)
    args = _base_args(
        tmp_path,
        output,
        _module_config(tmp_path),
        references,
        dry_run=True,
        mixer=str(fake_mixer),
        mixer_backend="native",
    )
    evidence = preflight_mixer_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )

    assert isinstance(evidence.resources, MixerPipelineResources)
    assert evidence.resources.backend == "native"
    assert {identity.path for identity in evidence.resources.file_identities} == {
        Path(references["bim"].replace("@", "1")).resolve(),
        Path(references["ld"].replace("@", "1")).resolve(),
        fake_mixer.resolve(),
        Path(sys.executable).resolve(),
    }

    captured = {}

    def fake_service(current_args, ctx, **kwargs):
        captured["input"] = current_args.mixer_input_file
        captured["resources"] = kwargs["pipeline_resources"]
        return {"analysis": "univariate"}

    monkeypatch.setattr(
        "postgwas.modules.mixer.service.run_mixer_direct",
        fake_service,
    )
    args._step_num = "04"
    ctx = RunContext(
        {"formatter": {"mixer": {"mixer_input": "formatted.sumstats.gz"}}},
        validations={"mixer": evidence},
    )

    assert run_mixer_runner(args, ctx) == {"analysis": "univariate"}
    assert captured == {
        "input": "formatted.sumstats.gz",
        "resources": evidence.resources,
    }
    assert args.output_directory == str(output)


def test_pipeline_mixer_rejects_reference_changed_after_preflight(tmp_path):
    references = _references(tmp_path)
    output = tmp_path / "out"
    args = _base_args(
        tmp_path,
        output,
        _module_config(tmp_path),
        references,
        dry_run=True,
        mixer_backend="docker",
    )
    evidence = preflight_mixer_pipeline(
        args,
        preflight_evidence=pipeline_input_vcf_evidence(),
    )
    Path(references["bim"].replace("@", "1")).write_text(
        "reference changed after validation\n",
        encoding="utf-8",
    )

    with pytest.raises(MixerError, match="changed after pipeline preflight"):
        run_mixer_direct(args, pipeline_resources=evidence.resources)

    assert not (output / "results").exists()
    assert "changed after pipeline preflight" in _canonical_log(output).read_text(
        encoding="utf-8",
    )


def test_univariate_dry_run_logs_exactly_fit1_and_test1(tmp_path):
    references = _references(tmp_path)
    output = tmp_path / "out"
    args = _base_args(
        tmp_path, output, _module_config(tmp_path), references,
        dry_run=True, threads=2, seed=17, mixer_backend="docker",
    )

    result = run_mixer_direct(args)

    univariate = result["univariate"]
    assert univariate["fit"].endswith("STUDY.fit.json")
    assert univariate["test"].endswith("STUDY.test.json")
    assert result["execution_backend"] == "docker"
    text = _canonical_log(output).read_text(encoding="utf-8")
    assert "ghcr.io/precimed/gsa-mixer:2.2.1" in text
    assert "fit1" in text and "test1" in text
    assert "--threads, 2" in text and "--seed, 17" in text
    assert "--chr2use, 1" in text and "status=COMPLETED" in text
    assert "trait2" not in text and "fit2" not in text and "test2" not in text
    resolved = yaml.safe_load(
        (output / "run_metadata" / "resolved_config.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert list(resolved["modules"]) == ["mixer"]
    assert set(resolved["resources"]) == {"executables", "containers"}
    assert set(resolved["resources"]["executables"]) == {
        "python", "mixer", "mixer_figures",
    }


def test_univariate_runner_writes_compact_results_and_official_figures(
    tmp_path, capsys,
):
    references = _references(tmp_path)
    fake_mixer = tmp_path / "fake_mixer.py"
    fake_figures = tmp_path / "fake_mixer_figures.py"
    _write_fake_mixer(fake_mixer)
    _write_fake_figures(fake_figures)
    output = tmp_path / "out"
    args = _base_args(
        tmp_path, output, _module_config(tmp_path), references,
        mixer=str(fake_mixer), mixer_figures=str(fake_figures),
    )

    result = run_mixer_direct(args)

    univariate = result["univariate"]
    fit = output / "results" / "raw" / "univariate" / "STUDY.fit.json"
    test = output / "results" / "raw" / "univariate" / "STUDY.test.json"
    assert json.loads(fit.read_text())["options"]["seed"] == 31
    assert json.loads(test.read_text())["options"]["seed"] == 31
    assert univariate["fit"] == str(fit) and univariate["test"] == str(test)
    assert result["execution_backend"] == "native"
    assert univariate["quality_status"] == "warning"
    summary = yaml.safe_load(Path(univariate["summary_yaml"]).read_text())
    assert summary["architecture"]["estimated_causal_variants"] == 20000.0
    assert summary["model_fit"]["mixture_vs_infinitesimal_delta_aic"] == 20.0
    assert summary["input_qc"]["variants_analysed"] == 1
    assert summary["diagnostics"]["qq_nonfinite_values"] == 1
    assert Path(univariate["summary_tsv"]).is_file()
    assert len(univariate["figures"]) == 4
    assert "--statistic point_estimate" in (output / "plots" / "STUDY.csv").read_text()
    screen = capsys.readouterr().out
    for milestone in (
        "Completed 1/4 · Validate the MiXeR input and resolved resources",
        "Completed 2/4 · Fit the single-trait MiXeR model",
        "Completed 3/4 · Evaluate the fitted single-trait MiXeR model",
        "Completed 4/4 · Validate and report the selected MiXeR results",
        "MiXeR fit1 execution progress",
        "MiXeR chromosome LD-loading progress",
        "Completed 1/1 · Load chromosome LD references",
        "MiXeR optimizer activity",
        "Completed 1/1 · Observed MiXeR cost-function evaluations",
    ):
        assert milestone in screen
    canonical_log = _canonical_log(output).read_text(encoding="utf-8")
    assert "mixer_optimizer_progress" in canonical_log
    assert "total_cost_evaluations=unknown_until_convergence" in canonical_log
    assert "mixer_native_progress_complete" in canonical_log


def test_univariate_native_failure_stops_progress_below_completion(
    tmp_path, capsys,
):
    references = _references(tmp_path)
    fake_mixer = tmp_path / "failing_mixer.py"
    _write_failing_mixer(fake_mixer)
    output = tmp_path / "out"
    args = _base_args(
        tmp_path,
        output,
        _module_config(tmp_path),
        references,
        mixer=str(fake_mixer),
    )

    with pytest.raises(MixerError, match="MiXeR fit1 failed with exit status 2"):
        run_mixer_direct(args)

    screen = capsys.readouterr().out
    assert "Failed 1/? · Observed MiXeR cost-function evaluations" in screen
    assert "Failed 3/3 · Optimize and validate the MiXeR model" in screen
    assert "Failed 2/4 · Fit the single-trait MiXeR model" in screen
    assert "Completed 2/4 · Fit the single-trait MiXeR model" not in screen
    assert "All 4 stages completed" not in screen


def test_gsa_runner_uses_official_three_stage_single_trait_contract(
    tmp_path, capsys,
):
    references = _references(tmp_path)
    go = _go_files(tmp_path)
    fake_mixer = tmp_path / "fake_mixer.py"
    _write_fake_mixer(fake_mixer)
    config = _module_config(
        tmp_path,
        analysis="gsa",
        extra=(
            "reporting:\n  generate_figures: false\n"
            "gsa:\n  top_results: 2\n  adam_epoch: [2, 3]\n"
            "  adam_step: [0.01, 0.001]\n"
        ),
    )
    output = tmp_path / "out"
    args = _base_args(
        tmp_path, output, config, references,
        mixer=str(fake_mixer),
        gsa_annotation_file_pattern=references["annotation"],
        gsa_baseline_go_file=go["baseline"],
        gsa_model_go_file=go["model"],
        gsa_test_go_file=go["test"],
    )

    result = run_mixer_direct(args)

    assert result["analysis"] == "gsa" and "univariate" not in result
    gsa = result["gsa"]
    summary = yaml.safe_load(Path(gsa["summary_yaml"]).read_text())
    assessment = summary["gene_set_assessment"]
    assert assessment["tested_identifiers"] == 3
    assert assessment["identifiers_with_positive_aic_evidence"] == 2
    assert summary["interpretation"]["multiple_testing_p_values_computed"] is False
    with Path(gsa["top_results_tsv"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["identifier"] for row in rows] == ["SET_C", "SET_A"]
    assert "p_value" not in rows[0] and "fdr" not in rows[0]
    assert all(Path(path).is_file() for path in (
        gsa["baseline_json"], gsa["full_json"], gsa["enrichment_results"],
    ))
    log_text = _canonical_log(output).read_text(encoding="utf-8")
    for token in (
        "split_sumstats", "--gsa-base", "--gsa-full", "--go-file-test",
        "--hardprune-maf, 0.05", "--hardprune-r2, 0.6",
        "--go-extend-bp, 10000", "--se-samples, 100",
        "--adam-epoch, 2, 3", "--adam-step, 0.01, 0.001",
    ):
        assert token in log_text
    assert "trait2" not in log_text and "fit2" not in log_text and "test2" not in log_text
    screen = capsys.readouterr().out
    for milestone in (
        "Completed 1/5 · Validate the MiXeR input and resolved resources",
        "Completed 2/5 · Split summary statistics by chromosome for GSA-MiXeR",
        "Completed 3/5 · Fit the GSA-MiXeR baseline model",
        "Completed 4/5 · Fit the GSA-MiXeR enrichment model",
        "Completed 5/5 · Validate and report the selected MiXeR results",
    ):
        assert milestone in screen


def test_gsa_precomputed_load_library_does_not_require_ld_files(tmp_path):
    references = _references(tmp_path, loadlib=True)
    go = _go_files(tmp_path)
    config = _module_config(tmp_path, analysis="gsa")
    output = tmp_path / "out"
    args = _base_args(
        tmp_path, output, config, references,
        dry_run=True, ld_file_pattern=None,
        gsa_annotation_file_pattern=references["annotation"],
        gsa_loadlib_file_pattern=references["loadlib"],
        gsa_baseline_go_file=go["baseline"],
        gsa_model_go_file=go["model"],
        gsa_test_go_file=go["test"],
    )

    run_mixer_direct(args)

    log_text = _canonical_log(output).read_text(encoding="utf-8")
    assert "--loadlib-file" in log_text
    assert "--ld-file" not in log_text


def test_gsa_missing_chromosome_resource_fails_before_external_execution(tmp_path):
    references = _references(tmp_path)
    go = _go_files(tmp_path)
    output = tmp_path / "out"
    args = _base_args(
        tmp_path, output, _module_config(tmp_path, analysis="gsa"), references,
        dry_run=True,
        gsa_annotation_file_pattern=str(tmp_path / "missing_chr@.annot.gz"),
        gsa_baseline_go_file=go["baseline"],
        gsa_model_go_file=go["model"],
        gsa_test_go_file=go["test"],
    )

    with pytest.raises(MixerError, match="Missing or empty GSA annotation"):
        run_mixer_direct(args)

    log_text = _canonical_log(output).read_text(encoding="utf-8")
    assert "FAILED" in log_text and "Missing or empty GSA annotation" in log_text
    assert "external_command" not in log_text


@pytest.mark.parametrize("invalid", ["fit_arguments: ['--trait2-file']\n", "replicates: 2\n"])
def test_configuration_rejects_multitrait_and_removed_legacy_options(tmp_path, invalid):
    config = _module_config(tmp_path, extra=invalid)
    with pytest.raises(ConfigurationError):
        load_module_configuration("mixer", config)


@pytest.mark.parametrize(
    "invalid",
    [
        "gsa:\n  evidence_aic_threshold: .inf\n",
        "gsa:\n  z_max: .inf\n",
        "gsa:\n  adam_epoch: [2]\n  adam_step: [.inf]\n",
    ],
)
def test_configuration_rejects_nonfinite_gsa_values(tmp_path, invalid):
    config = _module_config(tmp_path, extra=invalid)
    with pytest.raises(ConfigurationError):
        load_module_configuration("mixer", config)


def test_reference_failure_and_missing_input_are_always_logged(tmp_path):
    references = _references(tmp_path)
    output = tmp_path / "reference_failure"
    references["bim"] = str(tmp_path / "missing_chr@.bim")
    args = _base_args(tmp_path, output, _module_config(tmp_path), references)
    with pytest.raises(MixerError, match="Missing or empty BIM"):
        run_mixer_direct(args)
    assert "FAILED" in _canonical_log(output).read_text(encoding="utf-8")

    missing_output = tmp_path / "missing_input"
    with pytest.raises(MixerError, match="--mixer-input-file"):
        run_mixer_direct(Namespace(
            dataset_id="STUDY", output_directory=str(missing_output),
            dry_run=False, resume=False, overwrite=False,
        ))
    assert "FAILED" in _canonical_log(missing_output).read_text(encoding="utf-8")


def test_configured_nested_result_paths_are_preserved(tmp_path):
    references = _references(tmp_path)
    config = _module_config(
        tmp_path,
        extra=(
            "output_layout:\n"
            "  fit_prefix: 'fits/{dataset_id}.fit'\n"
            "  test_prefix: 'tests/{dataset_id}.test'\n"
        ),
    )
    output = tmp_path / "out"
    result = run_mixer_direct(_base_args(
        tmp_path, output, config, references, dry_run=True,
    ))
    assert result["univariate"]["fit"] == str(output / "fits" / "STUDY.fit.json")
    assert result["univariate"]["test"] == str(output / "tests" / "STUDY.test.json")
