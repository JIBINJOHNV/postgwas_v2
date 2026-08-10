"""Scientific and orchestration contracts for K-POPS and CALDERA."""

from __future__ import annotations

import argparse
from argparse import Namespace
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

from postgwas.config import load_configuration, load_module_configuration
from postgwas.modules.caldera.cli import build_parser as build_caldera_parser
from postgwas.modules.caldera.errors import CalderaError
from postgwas.modules.caldera.service import (
    convert_finemap_credible_sets,
    run_caldera_direct,
    validate_caldera_configuration,
)
from postgwas.modules.kpops.cli import build_parser as build_kpops_parser
from postgwas.modules.kpops.errors import KPopsError
from postgwas.modules.kpops.service import run_kpops_direct, validate_kpops_configuration
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY


def _kpops_resources(tmp_path: Path, *, kernel_bytes: int | None = None):
    genes = ["ENSG%03d" % number for number in range(1, 11)]
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "ENSGID\tNAME\tCHR\tTSS\n"
        + "".join(
            "%s\tGENE%d\t%d\t%d\n" % (gene, index, 1 + index % 2, 1000 + index)
            for index, gene in enumerate(genes)
        ),
        encoding="utf-8",
    )
    kernel_prefix = tmp_path / "kernel"
    Path(str(kernel_prefix) + ".genes").write_text("\n".join(genes) + "\n", encoding="utf-8")
    expected_bytes = len(genes) ** 2 * 4
    Path(str(kernel_prefix) + ".bin").write_bytes(b"\0" * (kernel_bytes or expected_bytes))
    magma_prefix = tmp_path / "magma"
    Path(str(magma_prefix) + ".genes.out").write_text(
        "GENE\tZSTAT\n"
        + "".join("%s\t%.3f\n" % (gene, index / 10) for index, gene in enumerate(genes)),
        encoding="utf-8",
    )
    Path(str(magma_prefix) + ".genes.raw").write_text(
        "header\nheader2\n"
        + "".join(
            "%s %d 0 0 10 5 0 100 0\n" % (gene, 1 if index < 5 else 2)
            for index, gene in enumerate(genes)
        ),
        encoding="utf-8",
    )
    script = tmp_path / "k-pops.py"
    script.write_text("# fixture\n", encoding="utf-8")
    return {
        "genes": genes, "annotation": annotation, "kernel_prefix": kernel_prefix,
        "magma_prefix": magma_prefix, "script": script,
    }


def _kpops_args(tmp_path: Path, resources, **overrides):
    values = {
        "dataset_id": "STUDY",
        "output_directory": str(tmp_path / "kpops_output"),
        "kpops_genome_build": "GRCh37",
        "kpops_script": str(resources["script"]),
        "kpops_gene_annotation_file": str(resources["annotation"]),
        "kernel_matrix_prefix": str(resources["kernel_prefix"]),
        "magma_association_prefix": str(resources["magma_prefix"]),
    }
    values.update(overrides)
    return Namespace(**values)


def _caldera_resources(tmp_path: Path):
    repository = tmp_path / "caldera"
    (repository / "data").mkdir(parents=True)
    (repository / "trained_models").mkdir()
    files = {
        "upstream": repository / "z_caldera.R",
        "coding": repository / "data" / "high_h2_coding_SNPs.tsv",
        "genes37": repository / "data" / "gene_locations_gencode_v47_grch37_all_protein_coding.tsv",
        "genes38": repository / "data" / "gene_locations_gencode_v47_grch38_all_protein_coding.tsv",
        "model": repository / "trained_models" / "caldera_model_no_covs.rds",
    }
    for path in files.values():
        path.write_text("fixture\n", encoding="utf-8")
    adapter = tmp_path / "run_caldera.R"
    adapter.write_text("# fixture\n", encoding="utf-8")
    pops = tmp_path / "pops.preds"
    pops.write_text("ENSGID\tPoPS_Score\nENSG001\t1.5\nENSG002\t-0.5\n", encoding="utf-8")
    credible_sets = tmp_path / "credible_sets.tsv"
    credible_sets.write_text(
        "locus\tchr\tbp\tpip\n"
        "locus1\t1\t100\t0.7\n"
        "locus1\t1\t110\t0.3\n",
        encoding="utf-8",
    )
    return {
        "repository": repository, "adapter": adapter, "pops": pops,
        "credible_sets": credible_sets,
    }


def _caldera_args(tmp_path: Path, resources, **overrides):
    values = {
        "dataset_id": "STUDY",
        "output_directory": str(tmp_path / "caldera_output"),
        "caldera_genome_build": "GRCh37",
        "caldera_repository": str(resources["repository"]),
        "caldera_adapter_script": str(resources["adapter"]),
        "pops_file": str(resources["pops"]),
        "credible_set_file": str(resources["credible_sets"]),
    }
    values.update(overrides)
    return Namespace(**values)


def test_clis_have_no_independent_defaults():
    for parser in (build_kpops_parser(), build_caldera_parser()):
        for action in parser._actions:
            if action.dest != "help":
                assert action.default == argparse.SUPPRESS
        assert vars(parser.parse_args([])) == {}


def test_kpops_data_paths_have_no_defaults_and_software_path_accepts_cli(tmp_path):
    defaults = load_configuration()
    assert defaults.modules.kpops.script_path == "k-pops.py"
    assert defaults.modules.kpops.gene_annotation_file is None
    assert defaults.modules.kpops.kernel_matrix_prefix is None

    resources = _kpops_resources(tmp_path)
    configuration, _ = validate_kpops_configuration(
        _kpops_args(tmp_path, resources)
    )
    assert configuration.modules.kpops.script_path == str(resources["script"])
    assert configuration.modules.kpops.gene_annotation_file == str(resources["annotation"])
    assert configuration.modules.kpops.kernel_matrix_prefix == str(
        resources["kernel_prefix"]
    )


def test_kpops_infers_installed_command_when_script_override_is_omitted(
    tmp_path, monkeypatch,
):
    resources = _kpops_resources(tmp_path)
    resources["script"].chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    args = _kpops_args(tmp_path, resources)
    del args.kpops_script

    configuration, observed = validate_kpops_configuration(args)

    assert configuration.modules.kpops.script_path == "k-pops.py"
    assert observed["script"] == resources["script"].resolve()


def test_registry_and_planner_connect_both_direct_and_pipeline_modes():
    assert REGISTRY.commands()["kpops"].cli_entrypoint.endswith("kpops.cli:main")
    assert REGISTRY.commands()["caldera"].cli_entrypoint.endswith("caldera.cli:main")
    assert build_pipeline_plan(["kpops"]).steps == ("formatter", "magma", "kpops")
    caldera = build_pipeline_plan(["caldera"]).steps
    assert caldera == (
        "annot_ldblock", "formatter", "ld_clump", "magma", "pops", "finemap", "caldera",
    )


def test_configuration_rejects_invalid_scientific_values(tmp_path):
    invalid_device = tmp_path / "kpops.yaml"
    invalid_device.write_text("device: gpu\n", encoding="utf-8")
    with pytest.raises(Exception, match="device"):
        load_module_configuration("kpops", invalid_device)
    invalid_coverage = tmp_path / "caldera.yaml"
    invalid_coverage.write_text("minimum_credible_set_coverage: 0.9\n", encoding="utf-8")
    with pytest.raises(Exception, match="minimum_credible_set_coverage"):
        load_module_configuration("caldera", invalid_coverage)


def test_kpops_validates_exact_kernel_dimensions_and_gene_overlap(tmp_path):
    resources = _kpops_resources(tmp_path)
    _, observed = validate_kpops_configuration(_kpops_args(tmp_path, resources))
    assert observed["kernel_gene_count"] == 10

    Path(str(resources["kernel_prefix"]) + ".bin").write_bytes(b"undersized")
    with pytest.raises(KPopsError, match="expected 400"):
        validate_kpops_configuration(_kpops_args(tmp_path, resources))


def test_kpops_rejects_unknown_anchor_genes(tmp_path):
    resources = _kpops_resources(tmp_path)
    with pytest.raises(KPopsError, match="anchor genes"):
        validate_kpops_configuration(
            _kpops_args(tmp_path, resources, anchor_genes=["NOT_A_GENE"])
        )


def test_kpops_allows_duplicate_symbols_but_rejects_ambiguous_name_anchor(tmp_path):
    resources = _kpops_resources(tmp_path)
    annotation = resources["annotation"]
    annotation.write_text(
        annotation.read_text(encoding="utf-8").replace("GENE1", "GENE0", 1),
        encoding="utf-8",
    )
    validate_kpops_configuration(_kpops_args(tmp_path, resources))
    with pytest.raises(KPopsError, match="map to multiple Ensembl IDs"):
        validate_kpops_configuration(
            _kpops_args(
                tmp_path, resources,
                anchor_gene_type="NAME", anchor_genes=["GENE0"],
            )
        )


def test_kpops_runs_in_staging_normalizes_gene_id_and_reports_findings(
    tmp_path, monkeypatch, capsys,
):
    resources = _kpops_resources(tmp_path)

    def fake_run(command, purpose, **kwargs):
        prefix = Path(command[command.index("--out_prefix") + 1])
        predictions = pd.DataFrame(
            {
                "NAME": ["GENE%d" % index for index in range(10)],
                "PoPS_Score": np.linspace(-1, 1, 10),
                "support_ENSGID": ["ENSG002,ENSG003"] + [None] * 9,
                "support_NAME": ["GENE1,GENE2"] + [None] * 9,
                "Anchor_Score": np.zeros(10),
            },
            index=resources["genes"],
        )
        predictions.to_csv(str(prefix) + ".preds", sep="\t")
        Path(str(prefix) + ".coefs").write_text("NAME\ttraining-1\nGENE1\t1\n", encoding="utf-8")
        return ""

    monkeypatch.setattr("postgwas.modules.kpops.service.run_checked_command", fake_run)
    result = run_kpops_direct(_kpops_args(tmp_path, resources))
    screen = capsys.readouterr().out
    published = pd.read_csv(result["kpops_file"], sep="\t")
    assert published.columns[0] == "ENSGID"
    assert "PoPS_Score" in published.columns
    assert published.loc[0, "support_ENSGID"] == "ENSG002,ENSG003"
    assert pd.isna(published.loc[1, "support_ENSGID"])
    assert len(published) == 10
    assert result["summary"]["genes_scored"] == 10
    assert result["summary"]["top_genes"][0]["gene_id"] == "ENSG010"
    assert "K-POPS gene-prioritisation summary" in screen
    assert "Genes with finite K-POPS scores" in screen
    assert "None. K-POPS scores are relative rankings, not p-values." in screen
    assert "Complete K-POPS results" in screen
    assert Path(result["completion_manifest"]).is_file()
    assert not (Path(_kpops_args(tmp_path, resources).output_directory) / ".partial").exists()

    resumed = run_kpops_direct(_kpops_args(tmp_path, resources))
    resumed_screen = capsys.readouterr().out
    assert resumed["summary"] == result["summary"]
    assert "K-POPS gene-prioritisation summary" in resumed_screen

    genes_out = Path(str(resources["magma_prefix"]) + ".genes.out")
    genes_out.write_text(
        genes_out.read_text(encoding="utf-8").replace("0.000", "0.001", 1),
        encoding="utf-8",
    )
    with pytest.raises(KPopsError, match="input magma_genes_out changed"):
        run_kpops_direct(_kpops_args(tmp_path, resources))


def test_kpops_summary_warns_about_incomplete_gene_and_chromosome_coverage(
    tmp_path, monkeypatch, capsys,
):
    resources = _kpops_resources(tmp_path)
    annotation = resources["annotation"]
    annotation.write_text(
        annotation.read_text(encoding="utf-8")
        .replace("ENSG009\tGENE8\t1\t1008", "ENSG009\tGENE8\t3\t1008")
        .replace("ENSG010\tGENE9\t2\t1009", "ENSG010\tGENE9\t3\t1009"),
        encoding="utf-8",
    )
    Path(str(resources["magma_prefix"]) + ".genes.out").write_text(
        "GENE\tZSTAT\n"
        + "".join(
            "%s\t%.3f\n" % (gene, index / 10)
            for index, gene in enumerate(resources["genes"][:8])
        ),
        encoding="utf-8",
    )
    Path(str(resources["magma_prefix"]) + ".genes.raw").write_text(
        "header\nheader2\n"
        + "".join(
            "%s %d 0 0 10 5 0 100 0\n"
            % (gene, 1 if index < 4 else 2)
            for index, gene in enumerate(resources["genes"][:8])
        ),
        encoding="utf-8",
    )
    run_config = tmp_path / "kpops_reporting.yaml"
    run_config.write_text(
        "minimum_gene_count: 2\n"
        "top_contributor_gene_count: 3\n"
        "reporting:\n"
        "  top_gene_count: 3\n"
        "  score_decimal_places: 3\n",
        encoding="utf-8",
    )

    def fake_run(command, purpose, **kwargs):
        prefix = Path(command[command.index("--out_prefix") + 1])
        predictions = pd.DataFrame(
            {
                "NAME": ["GENE%d" % index for index in range(10)],
                "PoPS_Score": list(np.linspace(-1, 1, 8)) + [None, None],
                "support_ENSGID": ["ENSG001"] * 8 + [None, None],
                "support_NAME": ["GENE0"] * 8 + [None, None],
                "Anchor_Score": np.zeros(10),
            },
            index=resources["genes"],
        )
        predictions.to_csv(str(prefix) + ".preds", sep="\t")
        Path(str(prefix) + ".coefs").write_text(
            "NAME\ttraining-1\nGENE1\t1\n", encoding="utf-8",
        )
        return ""

    monkeypatch.setattr(
        "postgwas.modules.kpops.service.run_checked_command", fake_run,
    )
    args = _kpops_args(tmp_path, resources)
    args.run_config = str(run_config)
    result = run_kpops_direct(args)

    screen = capsys.readouterr().out
    assert result["summary"]["genes_scored"] == 8
    assert len(result["summary"]["top_genes"]) == 3
    assert "COMPLETED WITH SCIENTIFIC WARNINGS" in screen
    assert "K-POPS scores cover 8 of 10 compatible kernel genes" in screen
    assert "No finite K-POPS scores were available on chromosome(s): 3." in screen


def test_caldera_rejects_missing_values_and_insufficient_pip_mass(tmp_path):
    resources = _caldera_resources(tmp_path)
    validate_caldera_configuration(_caldera_args(tmp_path, resources))
    resources["credible_sets"].write_text(
        "locus\tchr\tbp\tpip\n"
        "locus1\t1\t100\t0.5\n"
        "locus1\t1\t110\t0.44\n",
        encoding="utf-8",
    )
    with pytest.raises(CalderaError, match="below"):
        validate_caldera_configuration(_caldera_args(tmp_path, resources))
    resources["credible_sets"].write_text(
        "locus\tchr\tbp\tpip\n"
        "locus1\t1\t100\t0.7\n"
        "locus1\t1\t\t0.3\n",
        encoding="utf-8",
    )
    with pytest.raises(CalderaError, match="cannot be missing"):
        validate_caldera_configuration(_caldera_args(tmp_path, resources))


def test_caldera_accepts_exact_configured_credible_set_coverage(tmp_path):
    resources = _caldera_resources(tmp_path)
    resources["credible_sets"].write_text(
        "locus\tchr\tbp\tpip\n"
        "locus1\t1\t100\t0.5\n"
        "locus1\t1\t110\t0.45\n",
        encoding="utf-8",
    )
    _, observed = validate_caldera_configuration(_caldera_args(tmp_path, resources))
    assert observed["credible_set_metrics"]["upstream_retained_rows"] == 2


def test_caldera_infers_installed_repository_and_packaged_adapter(
    tmp_path, monkeypatch,
):
    resources = _caldera_resources(tmp_path)
    environment = tmp_path / "environment"
    python = environment / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("fixture\n", encoding="utf-8")
    python.chmod(0o755)
    installed = environment / "share" / "postgwas" / "caldera"
    shutil.copytree(resources["repository"], installed)
    monkeypatch.setenv("PATH", str(python.parent))
    args = _caldera_args(tmp_path, resources)
    del args.caldera_repository
    del args.caldera_adapter_script

    configuration, observed = validate_caldera_configuration(args)

    assert configuration.modules.caldera.repository_path is None
    assert configuration.modules.caldera.adapter_script_path is None
    assert observed["repository"] == installed.resolve()
    assert observed["adapter"] == (
        Path(__file__).parents[1]
        / "src" / "postgwas" / "modules" / "caldera" / "run_caldera.R"
    ).resolve()


def test_finemap_interchange_conversion_is_explicit_and_validated(tmp_path):
    resources = _caldera_resources(tmp_path)
    module = load_configuration(cli_overrides={
        "modules.caldera.repository_path": str(resources["repository"]),
        "modules.caldera.adapter_script_path": str(resources["adapter"]),
    }).modules.caldera
    interchange = tmp_path / "flames_input"
    interchange.mkdir()
    credible = interchange / "study_CS_L1.txt"
    credible.write_text(
        "index cred1 prob1\n"
        "1 1:100:A_G 0.7\n"
        "2 1:110:C_T 0.3\n",
        encoding="utf-8",
    )
    (interchange / "indexfile_rows.tsv").write_text(
        "Filename\tGenomicLocus\n"
        "study_CS_L1.txt\tchr1:1-200\n",
        encoding="utf-8",
    )
    destination = tmp_path / "converted.tsv"
    metrics = convert_finemap_credible_sets(interchange, destination, module)
    converted = pd.read_csv(destination, sep="\t")
    assert list(converted.columns) == ["locus", "chr", "bp", "pip"]
    assert converted["pip"].sum() == pytest.approx(1.0)
    assert metrics["loci"] == 1
    assert metrics["credible_set_files"] == 1


def test_caldera_direct_publishes_only_validated_results(tmp_path, monkeypatch):
    resources = _caldera_resources(tmp_path)

    def fake_run(command, purpose, **kwargs):
        output = Path(command[-2])
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            "locus\tlocus_pos\tgene\tcaldera\tmulti\tn_genes\tdist\tpops\tcoding\timpute_pops\tensgid\n"
            "locus1\t1_1_200\tGENE1\t0.6\t0.3\t2\t0\t1.5\t0\tFALSE\tENSG001\n"
            "locus1\t1_1_200\tGENE2\t0.4\t0.2\t2\t10\t-0.5\t0\tFALSE\tENSG002\n",
            encoding="utf-8",
        )
        return ""

    monkeypatch.setattr("postgwas.modules.caldera.service.run_checked_command", fake_run)
    result = run_caldera_direct(_caldera_args(tmp_path, resources))
    assert Path(result["caldera_file"]).is_file()
    assert Path(result["completion_manifest"]).is_file()


def test_docker_pins_upstream_commits_and_uses_main_environment():
    dockerfile = (Path(__file__).parents[1] / "Dockerfile").read_text(encoding="utf-8")
    environment = (Path(__file__).parents[1] / "environment.yml").read_text(encoding="utf-8")
    assert "8acd49ed8c96565b17c2997420f514f24b65e096" in dockerfile
    assert "81a8a0308741ae986660f711bbf6a7abbd3bca19" in dockerfile
    assert "install -m 0755 /opt/kpops/k-pops.py /opt/conda/envs/postgwas/bin/k-pops.py" in dockerfile
    assert "micromamba run -n postgwas python /opt/conda/envs/postgwas/bin/k-pops.py --help" in dockerfile
    assert "micromamba run -n postgwas Rscript" in dockerfile
    assert "mv /opt/caldera /opt/conda/envs/postgwas/share/postgwas/caldera" in dockerfile
    assert "pytorch>=2.1,<3" in environment


def test_resource_preparation_scripts_pin_and_verify_upstream_archives():
    root = Path(__file__).parents[1] / "tools" / "resource_preparation"
    preparation = (root / "prepare_kpops_caldera_resources.sh").read_text(
        encoding="utf-8",
    )
    smoke = (root / "test_kpops_caldera_resources.sh").read_text(
        encoding="utf-8",
    )
    for value in (
        "8acd49ed8c96565b17c2997420f514f24b65e096",
        "81a8a0308741ae986660f711bbf6a7abbd3bca19",
        "e162e01c59e084c3cb395c6f9171609ec535a2f6e769a213f042872f590745d3",
        "8ab5259671afe93767bca9db22e9088cb47d674b9cdf3733006fba09958a62b5",
    ):
        assert value in preparation
    assert "--resource-root" in preparation
    assert "--prepare-linear-kernel" in preparation
    assert "--python" in preparation
    assert 'install -m 0755 "${kpops_destination}/k-pops.py" "$kpops_command"' in preparation
    assert 'caldera_install="${python_environment}/share/postgwas/caldera"' in preparation
    assert "python -m postgwas" not in smoke
    assert "--kpops-script" not in smoke
    assert "--caldera-repository" not in smoke
    assert "--caldera-adapter-script" not in smoke
