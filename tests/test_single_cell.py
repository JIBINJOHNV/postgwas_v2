"""Scientific and orchestration regressions for single-cell integration."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from pathlib import Path
import shutil
import sys
import textwrap

import numpy as np
import pandas as pd
import pytest
import yaml

from postgwas.config import load_configuration
from postgwas.config.exporter import (
    render_module_configuration,
    render_pipeline_configuration,
)
from postgwas.config.models.modules.single_cell import SINGLE_CELL_TOOLS
from postgwas.modules.magmacovar.errors import MagmaCovarError
from postgwas.modules.single_cell.methods.magma_celltype.analysis import (
    normalize_magma_celltype_results,
    validate_magma_celltype_covariates,
)
from postgwas.modules.single_cell.methods.ldsc_celltype.analysis import (
    normalize_ldsc_celltype_results,
)
from postgwas.modules.single_cell.cli import build_parser
from postgwas.modules.single_cell.errors import SingleCellError
from postgwas.modules.single_cell.service import (
    preflight_single_cell_pipeline,
    resolve_single_cell_configuration,
    run_single_cell_direct,
)
from postgwas.modules.single_cell.methods.registry import METHOD_REGISTRY
from postgwas.modules.single_cell.methods.scdrs.runner import (
    validate_scdrs_covariates,
    validate_scdrs_h5ad,
)
from postgwas.pipeline.planner import build_pipeline_plan
from postgwas.pipeline.runners import run_formatter_runner, run_single_cell_runner


def _gene_results(path: Path, genes: int = 10) -> Path:
    rows = ["# VERSION = 110", "# COVAR = NSAMP MAC"]
    rows.extend(
        "%d 1 %d %d 20 5 1000 4.5 %.4f"
        % (gene, gene * 100, gene * 100 + 50, gene / 10)
        for gene in range(1, genes + 1)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _covariates(path: Path, *, include_average: bool = True) -> Path:
    header = "GENE Microglia Neuron"
    if include_average:
        header += " Average"
    rows = [header]
    for gene in range(1, 11):
        values = [str(gene), "%.3f" % (gene / 10), "%.3f" % (gene / 20)]
        if include_average:
            values.append("%.3f" % (gene / 15))
        rows.append(" ".join(values))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _gsa_results(path: Path) -> Path:
    path.write_text(
        "TOTAL_GENES = 10\n"
        "TEST_DIRECTION = greater (covar)\n"
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Microglia COVAR 10 0.2 0.1 0.05 0.01\n"
        "Neuron COVAR 10 0.1 0.05 0.04 0.04\n",
        encoding="utf-8",
    )
    return path


def _fake_magma(path: Path) -> Path:
    path.write_text(
        "#!%s\n" % sys.executable
        + "import pathlib\n"
        + "import sys\n"
        + "arguments = sys.argv[1:]\n"
        + "prefix = pathlib.Path(arguments[arguments.index('--out') + 1])\n"
        + "prefix.parent.mkdir(parents=True, exist_ok=True)\n"
        + "pathlib.Path(str(prefix) + '.log').write_text("
        + "' '.join(arguments) + '\\n', encoding='utf-8')\n"
        + "pathlib.Path(str(prefix) + '.gsa.out').write_text("
        + "'TOTAL_GENES = 10\\n'"
        + "'TEST_DIRECTION = greater (covar)\\n'"
        + "'VARIABLE TYPE NGENES BETA BETA_STD SE P\\n'"
        + "'Microglia COVAR 10 0.2 0.1 0.05 0.01\\n'"
        + "'Neuron COVAR 10 0.1 0.05 0.04 0.04\\n',"
        + " encoding='utf-8')\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _h5ad(
    path: Path,
    *,
    fractional: bool = False,
    gene_count: int = 300,
    empty_cells: int = 0,
) -> Path:
    anndata = pytest.importorskip("anndata")
    cells = ["cell_%03d" % index for index in range(60)]
    genes = ["GENE%03d" % index for index in range(gene_count)]
    value = 0.5 if fractional else 1
    matrix = np.full((len(cells), len(genes)), value, dtype=float)
    matrix[:empty_cells, :] = 0
    obs = pd.DataFrame(
        {
            "cell_type": ["Neuron"] * 30 + ["Microglia"] * 30,
            "state_score": np.linspace(0, 1, len(cells)),
            "unique_group": cells,
        },
        index=cells,
    )
    adata = anndata.AnnData(
        X=matrix,
        obs=obs,
        var=pd.DataFrame(index=genes),
        dtype=matrix.dtype,
    )
    adata.write_h5ad(path)
    return path


def _gene_sets(path: Path, *, genes: int = 51) -> Path:
    selected = ",".join("GENE%03d:1" % index for index in range(genes))
    path.write_text(
        "TRAIT\tGENESET\nSTUDY\t%s\n" % selected,
        encoding="utf-8",
    )
    return path


def _scdrs_covariates(
    path: Path, *, cells: int = 60, include_constant: bool = True,
) -> Path:
    columns = ["cell_id"]
    if include_constant:
        columns.append("const")
    columns.append("detected_genes")
    rows = ["\t".join(columns)]
    for index in range(cells):
        values = ["cell_%03d" % index]
        if include_constant:
            values.append("1")
        values.append(str(250 + index))
        rows.append("\t".join(values))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _magma_gene_statistics(path: Path, *, genes: int = 120) -> Path:
    rows = ["GENE CHR START STOP NSNPS NPARAM N ZSTAT P"]
    rows.extend(
        "%d 1 %d %d 5 2 10000 %.6f %.6g"
        % (
            index + 1,
            index * 100,
            index * 100 + 50,
            (genes - index) / 10,
            (index + 1) / (genes * 10),
        )
        for index in range(genes)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _identifier_map(path: Path, *, genes: int = 120) -> Path:
    rows = ["ENTREZID\tSYMBOL"]
    rows.extend(
        "%d\tGENE%03d" % (index + 1, index)
        for index in range(genes)
    )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _fake_scdrs(path: Path) -> Path:
    script = r"""
        #!{python}
        import csv
        import gzip
        import pathlib
        import sys

        args = sys.argv[1:]
        if args == ["scdrs", "__version__"]:
            print("1.0.3")
            raise SystemExit(0)
        command = args[0]
        if command == "munge-gs":
            source = pathlib.Path(args[args.index("--zscore-file") + 1])
            destination = pathlib.Path(args[args.index("--out-file") + 1])
            maximum = int(args[args.index("--n-max") + 1])
            with source.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                gene_column, trait = reader.fieldnames
                rows = list(reader)
            rows.sort(key=lambda row: float(row[trait]), reverse=True)
            genes = [row[gene_column] + ":" + row[trait] for row in rows[:maximum]]
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(
                "TRAIT\tGENESET\n" + trait + "\t" + ",".join(genes) + "\n",
                encoding="utf-8",
            )
            raise SystemExit(0)
        out = pathlib.Path(args[args.index("--out-folder") + 1])
        out.mkdir(parents=True, exist_ok=True)
        if command == "compute-score":
            gs = pathlib.Path(args[args.index("--gs-file") + 1])
            with gs.open(encoding="utf-8", newline="") as handle:
                traits = [
                    row["TRAIT"]
                    for row in csv.DictReader(handle, delimiter="\t")
                ]
            for trait in traits:
                for suffix in (".score.gz", ".full_score.gz"):
                    with gzip.open(out / (trait + suffix), "wt") as handle:
                        handle.write(
                            "index\traw_score\tnorm_score\tmc_pval\tpval\t"
                            "nlog10_pval\tzscore\n"
                        )
                        handle.write(
                            "cell_000\t1\t1\t0.5\t0.5\t0.301\t0\n"
                        )
        elif command == "perform-downstream":
            traits = [
                item.name[:-len(".full_score.gz")]
                for item in out.glob("*.full_score.gz")
            ]
            groups = (
                args[args.index("--group-analysis") + 1].split(",")
                if "--group-analysis" in args else []
            )
            for trait in traits:
                for group in groups:
                    (out / (trait + ".scdrs_group." + group)).write_text(
                        "group\tn_cell\tn_ctrl\tassoc_mcp\tassoc_mcz\t"
                        "hetero_mcp\thetero_mcz\n"
                        "Neuron\t30\t5\t0.1\t1\t0.5\t0\n",
                        encoding="utf-8",
                    )
                if "--corr-analysis" in args:
                    (out / (trait + ".scdrs_cell_corr")).write_text(
                        "variable\tn_cell\tcorr_mcp\tcorr_mcz\n"
                        "state_score\t60\t0.1\t1\n",
                        encoding="utf-8",
                    )
                if "--gene-analysis" in args:
                    (out / (trait + ".scdrs_gene")).write_text(
                        "index\tCORR\tRANK\nGENE000\t0.2\t0\n",
                        encoding="utf-8",
                    )
        """.format(python=sys.executable)
    path.write_text(textwrap.dedent(script).lstrip(), encoding="utf-8")
    path.chmod(0o755)
    return path


def _ldsc_reference_prefix(path: Path, *, include_m: bool = True) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    for chromosome in range(1, 23):
        Path("%s%d.l2.ldscore.gz" % (path, chromosome)).write_bytes(b"reference\n")
        if include_m:
            Path("%s%d.l2.M_5_50" % (path, chromosome)).write_text(
                "100\n", encoding="utf-8",
            )
    return str(path)


def _ldsc_inputs(tmp_path: Path) -> dict[str, str]:
    references = tmp_path / "references"
    baseline = _ldsc_reference_prefix(references / "baseline.")
    weights = _ldsc_reference_prefix(
        references / "weights.", include_m=False,
    )
    neuron = _ldsc_reference_prefix(references / "Neuron.")
    control = _ldsc_reference_prefix(references / "AllGenes.")
    ldcts = tmp_path / "brain.ldcts"
    ldcts.write_text(
        "Neuron\treferences/Neuron.,references/AllGenes.\n",
        encoding="utf-8",
    )
    sumstats = tmp_path / "study.sumstats.gz"
    with gzip.open(sumstats, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tZ\tN\nrs1\t2.5\t10000\n")
    return {
        "baseline": baseline,
        "weights": weights,
        "neuron": neuron,
        "control": control,
        "ldcts": str(ldcts),
        "sumstats": str(sumstats),
    }


def _fake_ldsc(path: Path) -> Path:
    script = r"""
        #!{python}
        import pathlib
        import sys

        args = sys.argv[1:]
        if not args:
            print("* Version 3.0.2")
            raise SystemExit(0)
        prefix = pathlib.Path(args[args.index("--out") + 1])
        ldcts = pathlib.Path(args[args.index("--ref-ld-chr-cts") + 1])
        labels = [line.split()[0] for line in ldcts.read_text().splitlines()]
        prefix.parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(str(prefix) + ".cell_type_results.txt").write_text(
            "Name\tCoefficient\tCoefficient_std_error\tCoefficient_P_value\n"
            + "".join(label + "\t0.2\t0.05\t0.01\n" for label in labels),
            encoding="utf-8",
        )
        pathlib.Path(str(prefix) + ".log").write_text(
            " ".join(args) + "\n", encoding="utf-8",
        )
        """.format(python=sys.executable)
    path.write_text(textwrap.dedent(script).lstrip(), encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_munge_sumstats(path: Path) -> Path:
    script = r"""
        #!{python}
        import gzip
        import pathlib
        import sys

        args = sys.argv[1:]
        prefix = pathlib.Path(args[args.index("--out") + 1])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(str(prefix) + ".sumstats.gz", "wt") as handle:
            handle.write("SNP\tZ\tN\nrs1\t2.5\t10000\n")
        pathlib.Path(str(prefix) + ".log").write_text(
            " ".join(args) + "\n", encoding="utf-8",
        )
        """.format(python=sys.executable)
    path.write_text(textwrap.dedent(script).lstrip(), encoding="utf-8")
    path.chmod(0o755)
    return path


def test_pipeline_plan_reuses_existing_magma_stage():
    plan = build_pipeline_plan(["single_cell"])
    assert plan.steps == ("formatter", "magma", "single_cell")
    assert "magmacovar" not in plan.steps


def test_ldsc_celltype_pipeline_plan_does_not_schedule_magma():
    plan = build_pipeline_plan(
        ["single_cell"],
        dependency_overrides={"single_cell": ("formatter",)},
    )

    assert plan.active_modules == ("formatter", "single_cell")
    assert plan.steps == ("formatter", "single_cell")


def test_packaged_configuration_declares_base_analysis_contract():
    module = load_configuration().modules.single_cell
    assert module.tools == ["magma_celltype"]
    assert module.magma_celltype.workflow == "base"
    assert module.magma_celltype.magmacovar_use_case == "tissue_specificity"
    assert module.magma_celltype.average_property == "Average"
    assert module.multiple_testing.methods == ["bonferroni", "fdr_bh"]
    assert module.scdrs.covariate_format.require_exact_cell_ids is True
    assert module.scdrs.covariate_format.minimum_cell_overlap_fraction == 0.75
    assert module.scdrs.covariate_format.constant_column == "const"


def test_scdrs_exact_mapping_requires_one_identifier_namespace():
    construction = (
        load_configuration()
        .modules.single_cell.scdrs.magma_gene_set.model_copy(deep=True)
    )

    with pytest.raises(ValueError, match="identical source and target"):
        construction.mapping_mode = "exact"


def test_scdrs_downstream_requires_normalized_control_scores():
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.return_control_raw_score = True
    method.downstream.group_analysis = ["cell_type"]

    with pytest.raises(ValueError, match="normalized control scores"):
        method.return_control_normalized_score = False


def test_scdrs_requires_identical_species_aliases():
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)

    with pytest.raises(ValueError, match="identical H5AD and gene-set species"):
        method.gene_set_species = "hsapiens"


def test_scdrs_annotation_names_reject_native_list_separator():
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)

    with pytest.raises(ValueError, match="must not contain commas"):
        method.downstream.group_analysis = ["cell,type"]


def test_tools_option_overrides_the_canonical_tool_selection():
    omitted = build_parser().parse_args([])
    args = build_parser().parse_args(["--tools", "magma_celltype"])
    configuration = resolve_single_cell_configuration(args)

    assert not hasattr(omitted, "tools")
    assert args.tools == ["magma_celltype"]
    assert configuration.modules.single_cell.tools == ["magma_celltype"]


def test_tools_option_accepts_scdrs():
    args = build_parser().parse_args(["--tools", "scdrs"])
    configuration = resolve_single_cell_configuration(args)

    assert args.tools == ["scdrs"]
    assert configuration.modules.single_cell.tools == ["scdrs"]


def test_tools_option_accepts_ldsc_celltype():
    args = build_parser().parse_args(["--tools", "ldsc_celltype"])
    configuration = resolve_single_cell_configuration(args)

    assert args.tools == ["ldsc_celltype"]
    assert configuration.modules.single_cell.tools == ["ldsc_celltype"]


def test_method_registry_matches_configuration_and_exposes_one_adapter_per_tool():
    assert tuple(METHOD_REGISTRY) == tuple(SINGLE_CELL_TOOLS)
    for name, method in METHOD_REGISTRY.items():
        assert method.name == name
        assert method.resource_paths
        assert callable(method.prepare_pipeline_args)
        assert callable(method.preflight_pipeline)
        assert callable(method.preflight_direct)
        assert callable(method.run)


def test_single_cell_help_discloses_scdrs_internal_seed():
    assert "fixed internal seed 0" in build_parser().format_help()


def test_pipeline_configuration_export_selects_magma_formatter_contract():
    document = yaml.safe_load(
        render_pipeline_configuration(["single_cell"], style="values")
    )

    assert list(document["modules"]) == [
        "formatting", "magma", "magmacovar", "single_cell",
    ]
    assert document["pipeline"]["modules"] == [
        "formatting", "magma", "single_cell",
    ]
    assert document["modules"]["formatting"]["formats"] == ["magma"]
    assert document["modules"]["magmacovar"]["enabled"] is False
    assert document["modules"]["single_cell"]["tools"] == ["magma_celltype"]


def test_ldsc_celltype_pipeline_export_selects_only_ldsc_formatter_contract(
    tmp_path,
):
    config = tmp_path / "pipeline.yaml"
    config.write_text(
        "modules:\n"
        "  single_cell:\n"
        "    tools: [ldsc_celltype]\n"
        "  ldsc:\n"
        "    minimum_info: 0.8\n",
        encoding="utf-8",
    )
    document = yaml.safe_load(render_pipeline_configuration(
        ["single_cell"], config_file=config, style="values",
    ))

    assert list(document["modules"]) == ["formatting", "ldsc", "single_cell"]
    assert document["pipeline"]["modules"] == ["formatting", "single_cell"]
    assert document["modules"]["formatting"]["formats"] == ["ldsc"]
    assert document["modules"]["ldsc"]["minimum_info"] == 0.8
    assert document["modules"]["single_cell"]["tools"] == ["ldsc_celltype"]


def test_covariate_validation_requires_average_expression(tmp_path):
    configuration = load_configuration()
    covariates = _covariates(
        tmp_path / "without_average.tsv", include_average=False,
    )
    with pytest.raises(SingleCellError, match="average expression property"):
        validate_magma_celltype_covariates(
            covariates,
            average_property="Average",
            magmacovar_config=configuration.modules.magmacovar,
        )


def test_normalization_corrects_across_the_complete_cell_type_family(tmp_path):
    configuration = load_configuration()
    output = tmp_path / "cell_types.tsv"
    metrics = normalize_magma_celltype_results(
        _gsa_results(tmp_path / "study.gsa.out"),
        _covariates(tmp_path / "covariates.tsv"),
        output,
        dataset_id="study",
        single_cell_config=configuration.modules.single_cell,
        magmacovar_config=configuration.modules.magmacovar,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["cell_type"] for row in rows] == ["Microglia", "Neuron"]
    assert {row["dataset_id"] for row in rows} == {"study"}
    assert [float(row["p_bonferroni"]) for row in rows] == [0.02, 0.08]
    assert [float(row["p_fdr_bh"]) for row in rows] == [0.02, 0.04]
    assert metrics["tested_cell_types"] == 2
    assert metrics["adjusted_significant"] == {"bonferroni": 1, "fdr_bh": 2}


def test_normalization_rejects_incomplete_magma_hypothesis_family(tmp_path):
    configuration = load_configuration()
    incomplete = tmp_path / "incomplete.gsa.out"
    incomplete.write_text(
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Microglia COVAR 10 0.2 0.1 0.05 0.01\n",
        encoding="utf-8",
    )
    with pytest.raises(SingleCellError, match="missing cell types: Neuron"):
        normalize_magma_celltype_results(
            incomplete,
            _covariates(tmp_path / "covariates.tsv"),
            tmp_path / "output.tsv",
            dataset_id="study",
            single_cell_config=configuration.modules.single_cell,
            magmacovar_config=configuration.modules.magmacovar,
        )


def test_normalization_rejects_duplicate_cell_type_results(tmp_path):
    configuration = load_configuration()
    duplicate = tmp_path / "duplicate.gsa.out"
    duplicate.write_text(
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Microglia COVAR 10 0.2 0.1 0.05 0.01\n"
        "Microglia COVAR 10 0.3 0.2 0.05 0.02\n"
        "Neuron COVAR 10 0.1 0.05 0.04 0.04\n",
        encoding="utf-8",
    )

    with pytest.raises(MagmaCovarError, match="duplicate COVAR variables"):
        normalize_magma_celltype_results(
            duplicate,
            _covariates(tmp_path / "covariates.tsv"),
            tmp_path / "output.tsv",
            dataset_id="study",
            single_cell_config=configuration.modules.single_cell,
            magmacovar_config=configuration.modules.magmacovar,
        )


def test_ldsc_celltype_normalization_preserves_native_one_sided_test(tmp_path):
    native = tmp_path / "study.cell_type_results.txt"
    native.write_text(
        "Name\tCoefficient\tCoefficient_std_error\tCoefficient_P_value\n"
        "Neuron\t0.2\t0.05\t0.01\n"
        "Microglia\t0.1\t0.04\t0.04\n",
        encoding="utf-8",
    )
    output = tmp_path / "normalized.tsv"
    metrics = normalize_ldsc_celltype_results(
        native,
        output,
        expected_cell_types=("Neuron", "Microglia"),
        dataset_id="study",
        single_cell_config=load_configuration().modules.single_cell,
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [float(row["beta"]) for row in rows] == [0.2, 0.1]
    assert [float(row["p_bonferroni"]) for row in rows] == [0.02, 0.08]
    assert rows[0]["n_genes"] == "NA"
    assert metrics["native_p_value_alternative"] == "coefficient_greater_than_zero"


def test_ldsc_celltype_direct_runs_validated_h2_cts_and_resumes(
    tmp_path, monkeypatch,
):
    inputs = _ldsc_inputs(tmp_path)
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        ldsc_celltype_sumstats_file=inputs["sumstats"],
        ldsc_celltype_ldcts_file=inputs["ldcts"],
        ldsc_celltype_baseline_prefix=[inputs["baseline"]],
        ldsc_celltype_weights_prefix=inputs["weights"],
        ldsc=str(_fake_ldsc(tmp_path / "ldsc.py")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )

    first = run_single_cell_direct(args)

    normalized = first.artifacts["ldsc_celltype_normalized_results"].path
    assert normalized.is_file()
    assert first.metrics["tested_cell_types"] == 1
    assert first.metrics["reference_files_validated"] == 154
    native_log = first.artifacts["ldsc_celltype_native_log"].path.read_text(
        encoding="utf-8"
    )
    assert "--h2-cts" in native_log
    assert "--ref-ld-chr-cts" in native_log
    qc = yaml.safe_load(
        first.artifacts["ldsc_celltype_qc_report"].path.read_text(
            encoding="utf-8"
        )
    )
    assert qc["scientific_settings"]["p_value_alternative"] == (
        "coefficient_greater_than_zero"
    )

    original = __import__(
        "postgwas.modules.single_cell.methods.ldsc_celltype.runner",
        fromlist=["run_checked_command"],
    ).run_checked_command

    def version_only(arguments, purpose, **kwargs):
        if purpose == "LDSC version probe":
            return original(arguments, purpose, **kwargs)
        pytest.fail("resume reran LDSC")

    monkeypatch.setattr(
        "postgwas.modules.single_cell.methods.ldsc_celltype.runner."
        "run_checked_command",
        version_only,
    )
    resumed = run_single_cell_direct(args)
    assert resumed.metrics["resumed"] is True
    assert resumed.artifacts["ldsc_celltype_normalized_results"].path == normalized


def test_ldsc_celltype_formatter_source_munges_before_h2_cts(tmp_path):
    inputs = _ldsc_inputs(tmp_path)
    raw_sumstats = tmp_path / "formatter_ldsc.tsv"
    raw_sumstats.write_text(
        "SNP\tA1\tA2\tZ\tP\tN\tFRQ\tINFO\n"
        "rs1\tG\tA\t2.5\t0.01\t10000\t0.2\t0.99\n",
        encoding="utf-8",
    )
    merge = tmp_path / "w_hm3.snplist"
    merge.write_text("SNP\nrs1\n", encoding="utf-8")
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        ldsc_celltype_sumstats_source="formatter",
        ldsc_celltype_sumstats_file=str(raw_sumstats),
        ldsc_celltype_ldcts_file=inputs["ldcts"],
        ldsc_celltype_baseline_prefix=[inputs["baseline"]],
        ldsc_celltype_weights_prefix=inputs["weights"],
        ldsc_celltype_merge_alleles_file=str(merge),
        ldsc=str(_fake_ldsc(tmp_path / "ldsc.py")),
        munge_sumstats=str(_fake_munge_sumstats(tmp_path / "munge_sumstats.py")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    result = run_single_cell_direct(args)

    assert result.metrics["sumstats_source"] == "formatter"
    assert result.artifacts["ldsc_celltype_munged_sumstats"].path.is_file()
    munge_log = result.artifacts["ldsc_celltype_munge_log"].path.read_text(
        encoding="utf-8"
    )
    assert "--info-min 0.9" in munge_log
    assert "--maf-min 0.01" in munge_log


def test_ldsc_celltype_requires_test_and_control_prefixes(tmp_path):
    inputs = _ldsc_inputs(tmp_path)
    Path(inputs["ldcts"]).write_text(
        "Neuron\treferences/Neuron.\n", encoding="utf-8",
    )
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        ldsc_celltype_sumstats_file=inputs["sumstats"],
        ldsc_celltype_ldcts_file=inputs["ldcts"],
        ldsc_celltype_baseline_prefix=[inputs["baseline"]],
        ldsc_celltype_weights_prefix=inputs["weights"],
        ldsc=str(_fake_ldsc(tmp_path / "ldsc.py")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="tested annotation followed"):
        run_single_cell_direct(args)


def test_ldsc_celltype_requires_complete_autosomal_reference_panel(tmp_path):
    inputs = _ldsc_inputs(tmp_path)
    Path(inputs["baseline"] + "22.l2.M_5_50").unlink()
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        ldsc_celltype_sumstats_file=inputs["sumstats"],
        ldsc_celltype_ldcts_file=inputs["ldcts"],
        ldsc_celltype_baseline_prefix=[inputs["baseline"]],
        ldsc_celltype_weights_prefix=inputs["weights"],
        ldsc=str(_fake_ldsc(tmp_path / "ldsc.py")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="22.l2.M_5_50"):
        run_single_cell_direct(args)


def test_direct_service_reuses_magmacovar_with_fuma_base_model(tmp_path):
    magma = _fake_magma(tmp_path / "magma")
    args = argparse.Namespace(
        magma_gene_results_file=str(_gene_results(tmp_path / "study.genes.raw")),
        single_cell_covariates=str(_covariates(tmp_path / "covariates.tsv")),
        magma=str(magma),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )

    result = run_single_cell_direct(args)

    normalized = result.artifacts["magma_celltype_results"].path
    assert normalized.is_file()
    engine_log = (
        tmp_path
        / "results"
        / "engines"
        / "study_magma_celltype"
        / "study.log"
    ).read_text(encoding="utf-8")
    assert "condition-hide=Average" in engine_log
    assert "direction-covar=greater" in engine_log
    assert (
        tmp_path
        / "results"
        / "run_metadata"
        / "study_single_cell_completion.yaml"
    ).is_file()


def test_direct_service_runs_registered_methods_in_selected_order(tmp_path):
    args = argparse.Namespace(
        tools=["magma_celltype", "scdrs"],
        magma_gene_results_file=str(
            _gene_results(tmp_path / "study.genes.raw")
        ),
        single_cell_covariates=str(
            _covariates(tmp_path / "covariates.tsv")
        ),
        magma=str(_fake_magma(tmp_path / "magma")),
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    result = run_single_cell_direct(args)

    assert list(result.metrics["tools"]) == ["magma_celltype", "scdrs"]
    assert result.artifacts["magma_celltype_results"].path.is_file()
    assert result.artifacts["scdrs_score_1"].path.is_file()
    assert result.artifacts["scdrs_full_score_1"].path.is_file()


def test_direct_service_accepts_exported_module_only_configuration(tmp_path):
    module_config = tmp_path / "single_cell.yaml"
    module_config.write_text(
        render_module_configuration("single_cell", style="values"),
        encoding="utf-8",
    )
    args = argparse.Namespace(
        magma_gene_results_file=str(_gene_results(tmp_path / "study.genes.raw")),
        single_cell_covariates=str(_covariates(tmp_path / "covariates.tsv")),
        magma=str(_fake_magma(tmp_path / "magma")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        run_config=str(module_config),
        resume=True,
    )

    result = run_single_cell_direct(args)

    assert result.artifacts["magma_celltype_results"].path.is_file()


def test_direct_service_resumes_delegated_and_normalized_results(
    tmp_path, monkeypatch,
):
    args = argparse.Namespace(
        magma_gene_results_file=str(_gene_results(tmp_path / "study.genes.raw")),
        single_cell_covariates=str(_covariates(tmp_path / "covariates.tsv")),
        magma=str(_fake_magma(tmp_path / "magma")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )
    first = run_single_cell_direct(args)

    monkeypatch.setattr(
        "postgwas.modules.magmacovar.service.run_magma_covariates",
        lambda **kwargs: pytest.fail("resume reran MAGMAcovar"),
    )
    monkeypatch.setattr(
        "postgwas.modules.single_cell.methods.magma_celltype.method."
        "normalize_magma_celltype_results",
        lambda *args, **kwargs: pytest.fail("resume renormalized cell types"),
    )
    args.scdrs_h5ad_file = str(tmp_path / "irrelevant_scdrs_atlas.h5ad")
    resumed = run_single_cell_direct(args)

    assert resumed.metrics["resumed"] is True
    assert (
        resumed.artifacts["magma_celltype_results"].path
        == first.artifacts["magma_celltype_results"].path
    )


def test_scdrs_direct_runs_native_scores_and_downstream_outputs(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        scdrs_group_analysis=["cell_type"],
        scdrs_correlation_analysis=["state_score"],
        scdrs_gene_analysis=True,
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )

    result = run_single_cell_direct(args)

    assert result.metrics["traits"] == 1
    assert result.metrics["software_version"] == "1.0.3"
    assert result.artifacts["scdrs_score_1"].path.is_file()
    assert result.artifacts["scdrs_full_score_1"].path.is_file()
    assert result.artifacts["scdrs_group_1_1"].path.is_file()
    assert result.artifacts["scdrs_correlation_1"].path.is_file()
    assert result.artifacts["scdrs_gene_1"].path.is_file()
    assert result.artifacts["scdrs_qc_report"].path.is_file()
    qc = yaml.safe_load(
        result.artifacts["scdrs_qc_report"].path.read_text(encoding="utf-8")
    )
    assert qc["metrics"]["upstream_internal_random_seed"] == 0
    assert qc["metrics"]["postgwas_execution_seed_applied"] is False
    annotations = qc["h5ad_validation"]["annotation_validation"]
    group_counts = annotations["analysis_group_counts"]["cell_type"]
    assert {row["value"]: row["cells"] for row in group_counts} == {
        "Neuron": 30,
        "Microglia": 30,
    }
    assert annotations["input_group_counts"] == annotations[
        "analysis_group_counts"
    ]
    assert (
        tmp_path
        / "results"
        / "run_metadata"
        / "study_scdrs_completion.yaml"
    ).is_file()


def test_scdrs_direct_resumes_without_rerunning_native_cli(tmp_path, monkeypatch):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )
    first = run_single_cell_direct(args)
    original = __import__(
        "postgwas.modules.single_cell.methods.scdrs.runner",
        fromlist=["run_checked_command"],
    ).run_checked_command

    def version_only(arguments, purpose, **kwargs):
        if "version probe" in purpose:
            return original(arguments, purpose, **kwargs)
        pytest.fail("resume reran native scDRS")

    monkeypatch.setattr(
        "postgwas.modules.single_cell.methods.scdrs.runner.run_checked_command",
        version_only,
    )
    resumed = run_single_cell_direct(args)

    assert resumed.metrics["resumed"] is True
    assert (
        resumed.artifacts["scdrs_score_1"].path
        == first.artifacts["scdrs_score_1"].path
    )


def test_scdrs_restarts_when_every_declared_output_is_missing(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )
    first = run_single_cell_direct(args)
    manifest = (
        tmp_path
        / "results"
        / "run_metadata"
        / "study_scdrs_completion.yaml"
    )
    recorded = yaml.safe_load(manifest.read_text(encoding="utf-8"))
    for fingerprint in recorded["outputs"].values():
        Path(fingerprint["path"]).unlink()

    restarted = run_single_cell_direct(args)

    assert restarted.metrics["resumed"] is False
    assert restarted.artifacts["scdrs_score_1"].path.is_file()
    log = tmp_path / "results" / "logs" / "study_single_cell.log"
    assert "reason=incomplete_outputs" in log.read_text(
        encoding="utf-8"
    )


def test_scdrs_qc_reports_group_counts_after_cell_filtering(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(
            _h5ad(tmp_path / "atlas.h5ad", empty_cells=10)
        ),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        scdrs_group_analysis=["cell_type"],
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    result = run_single_cell_direct(args)
    qc = yaml.safe_load(
        result.artifacts["scdrs_qc_report"].path.read_text(encoding="utf-8")
    )
    annotations = qc["h5ad_validation"]["annotation_validation"]
    input_counts = {
        row["value"]: row["cells"]
        for row in annotations["input_group_counts"]["cell_type"]
    }
    analysis_counts = {
        row["value"]: row["cells"]
        for row in annotations["analysis_group_counts"]["cell_type"]
    }
    assert input_counts == {"Neuron": 30, "Microglia": 30}
    assert analysis_counts == {"Neuron": 20, "Microglia": 30}


def test_scdrs_rejects_overgranular_proportion_adjustment(tmp_path):
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.adjust_proportion_column = "unique_group"

    with pytest.raises(SingleCellError, match="fewer than one group per ten"):
        validate_scdrs_h5ad(_h5ad(tmp_path / "atlas.h5ad"), method)


def test_scdrs_raw_count_declaration_rejects_fractional_matrix(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(
            _h5ad(tmp_path / "fractional.h5ad", fractional=True)
        ),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "study.gs")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="non-integer values"):
        run_single_cell_direct(args)


def test_scdrs_covariates_validate_constant_and_exact_cell_alignment(tmp_path):
    method = load_configuration().modules.single_cell.scdrs
    summary = validate_scdrs_covariates(
        _scdrs_covariates(tmp_path / "cells.cov"),
        method,
        cell_names=tuple("cell_%03d" % index for index in range(60)),
    )

    assert summary["overlap_fraction"] == 1
    assert summary["constant_column"] == "const"
    assert summary["maximum_constant_deviation"] == 0


def test_scdrs_covariates_require_configured_constant_column(tmp_path):
    method = load_configuration().modules.single_cell.scdrs
    with pytest.raises(SingleCellError, match="constant column 'const'"):
        validate_scdrs_covariates(
            _scdrs_covariates(
                tmp_path / "no_constant.cov", include_constant=False,
            ),
            method,
            cell_names=tuple("cell_%03d" % index for index in range(60)),
        )


def test_scdrs_covariates_require_more_than_configured_overlap(tmp_path):
    method = load_configuration().modules.single_cell.scdrs.model_copy(deep=True)
    method.covariate_format.require_exact_cell_ids = False
    with pytest.raises(SingleCellError, match="overlap 0.75"):
        validate_scdrs_covariates(
            _scdrs_covariates(tmp_path / "partial.cov", cells=45),
            method,
            cell_names=tuple("cell_%03d" % index for index in range(60)),
        )


def test_scdrs_rejects_gene_sets_below_effective_gene_minimum(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(
            _gene_sets(tmp_path / "too_small.gs", genes=50)
        ),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="configured minimum is 51"):
        run_single_cell_direct(args)


def test_scdrs_rejects_mixed_weighted_and_unweighted_gene_set(tmp_path):
    mixed = _gene_sets(tmp_path / "mixed.gs")
    mixed.write_text(
        mixed.read_text(encoding="utf-8").replace("GENE050:1", "GENE050"),
        encoding="utf-8",
    )
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(_h5ad(tmp_path / "atlas.h5ad")),
        scdrs_gene_set_file=str(mixed),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="mixes weighted and unweighted"):
        run_single_cell_direct(args)


def test_scdrs_constructs_pipeline_gene_set_from_magma_with_crosswalk(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_gene_set_source="magma",
        scdrs_h5ad_file=str(
            _h5ad(tmp_path / "atlas.h5ad", gene_count=1000)
        ),
        scdrs_magma_gene_results_file=str(
            _magma_gene_statistics(tmp_path / "study.genes.out")
        ),
        scdrs_gene_id_map=str(_identifier_map(tmp_path / "entrez_symbol.tsv")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
        resume=True,
    )

    result = run_single_cell_direct(args)

    assert result.metrics["gene_set_source"] == "magma"
    generated = result.artifacts["scdrs_generated_gene_set"].path
    statistics = result.artifacts["scdrs_generated_gene_statistics"].path
    mapping_report = result.artifacts["scdrs_gene_mapping_report"].path
    assert generated.is_file()
    assert statistics.is_file()
    assert mapping_report.is_file()
    assert generated.read_text(encoding="utf-8").startswith(
        "TRAIT\tGENESET\nstudy\tGENE000:"
    )
    report = yaml.safe_load(mapping_report.read_text(encoding="utf-8"))
    assert report["input_magma_genes"] == 120
    assert report["mapped_magma_genes"] == 120
    assert report["unmapped_magma_genes"] == 0
    assert report["source_identifier_type"] == "entrez"
    assert report["target_identifier_type"] == "symbol"


def test_scdrs_magma_crosswalk_rejects_ambiguous_target_ids(tmp_path):
    crosswalk = _identifier_map(tmp_path / "ambiguous.tsv")
    with crosswalk.open("a", encoding="utf-8") as handle:
        handle.write("121\tGENE000\n")
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_gene_set_source="magma",
        scdrs_h5ad_file=str(
            _h5ad(tmp_path / "atlas.h5ad", gene_count=1000)
        ),
        scdrs_magma_gene_results_file=str(
            _magma_gene_statistics(tmp_path / "study.genes.out")
        ),
        scdrs_gene_id_map=str(crosswalk),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    with pytest.raises(SingleCellError, match="not one-to-one"):
        run_single_cell_direct(args)


def test_scdrs_pipeline_preflight_defers_only_the_magma_result(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(
            _h5ad(tmp_path / "atlas.h5ad", gene_count=1000)
        ),
        scdrs_gene_id_map=str(_identifier_map(tmp_path / "entrez_symbol.tsv")),
        scdrs=str(_fake_scdrs(tmp_path / "scdrs")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    preflight_single_cell_pipeline(args)

    assert args.scdrs_gene_set_source == "magma"
    assert not hasattr(args, "scdrs_magma_gene_results_file")


def test_ldsc_celltype_pipeline_preflight_defers_only_formatter_sumstats(
    tmp_path,
):
    inputs = _ldsc_inputs(tmp_path)
    merge = tmp_path / "w_hm3.snplist"
    merge.write_text("SNP\nrs1\n", encoding="utf-8")
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        ldsc_celltype_ldcts_file=inputs["ldcts"],
        ldsc_celltype_baseline_prefix=[inputs["baseline"]],
        ldsc_celltype_weights_prefix=inputs["weights"],
        ldsc_celltype_merge_alleles_file=str(merge),
        ldsc=str(_fake_ldsc(tmp_path / "ldsc.py")),
        munge_sumstats=str(_fake_munge_sumstats(tmp_path / "munge_sumstats.py")),
        dataset_id="study",
        output_directory=str(tmp_path / "results"),
    )

    preflight_single_cell_pipeline(args)

    assert args.ldsc_celltype_sumstats_source == "formatter"
    assert not hasattr(args, "ldsc_celltype_sumstats_file")


def test_scdrs_pipeline_runner_injects_validated_magma_gene_output(
    tmp_path, monkeypatch,
):
    genes_out = _magma_gene_statistics(tmp_path / "study.genes.out")
    genes_raw = _gene_results(tmp_path / "study.genes.raw")
    args = argparse.Namespace(
        tools=["scdrs"],
        output_directory=str(tmp_path / "results"),
        dataset_id="study",
        _step_num=3,
    )
    context = {
        "magma": {
            "primary_mapping": "positional",
            "mapping_analyses": {
                "positional": {
                    "result_statistic_type": "calibrated_gene_p_value",
                    "gene_id_type": "entrez",
                },
            },
            "magma_genes_raw": str(genes_raw),
            "magma_genes_out": str(genes_out),
        },
    }
    observed = {}

    def fake_run(runtime_args, runtime_context):
        observed["source"] = runtime_args.scdrs_gene_set_source
        observed["genes"] = runtime_args.scdrs_magma_gene_results_file
        observed["context"] = runtime_context
        return "completed"

    monkeypatch.setattr(
        "postgwas.modules.single_cell.service.run_single_cell_direct",
        fake_run,
    )

    result = run_single_cell_runner(args, context)

    assert result == "completed"
    assert observed["source"] == "magma"
    assert observed["genes"] == str(genes_out)
    assert observed["context"] is context
    assert args.output_directory == str(tmp_path / "results")


def test_scdrs_pipeline_runner_rejects_magma_identifier_mismatch(tmp_path):
    args = argparse.Namespace(
        tools=["scdrs"],
        output_directory=str(tmp_path / "results"),
        dataset_id="study",
        _step_num=3,
    )
    context = {
        "magma": {
            "primary_mapping": "symbol_mapping",
            "mapping_analyses": {
                "symbol_mapping": {
                    "result_statistic_type": "calibrated_gene_p_value",
                    "gene_id_type": "symbol",
                },
            },
            "magma_genes_raw": "unused.genes.raw",
            "magma_genes_out": "unused.genes.out",
        },
    }

    with pytest.raises(ValueError, match="produced 'symbol'"):
        run_single_cell_runner(args, context)


def test_ldsc_celltype_pipeline_formatter_selects_ldsc_without_magma_reference(
    tmp_path, monkeypatch,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf\n")
    captured = {}

    def fake_formatter(args, context):
        captured["formats"] = list(args.format)
        captured["variant_id_types"] = dict(args.variant_id_types)
        captured["context"] = context
        return {"ldsc": {"ldsc_file": "study_ldsc.tsv"}}

    monkeypatch.setattr(
        "postgwas.modules.formatting.service.run_formatter_direct",
        fake_formatter,
    )
    args = argparse.Namespace(
        modules=["formatter", "single_cell"],
        tools=["ldsc_celltype"],
        vcf=str(vcf),
        bcftools=sys.executable,
        output_directory=str(tmp_path / "results"),
        dataset_id="study",
        _step_num=1,
    )
    context = {}

    result = run_formatter_runner(args, context)

    assert result == {"ldsc": {"ldsc_file": "study_ldsc.tsv"}}
    assert captured["formats"] == ["ldsc"]
    assert captured["variant_id_types"] == {"ldsc": "rsid"}
    assert captured["context"] is context
    assert args.output_directory == str(tmp_path / "results")


def test_ldsc_celltype_pipeline_rejects_non_rsid_formatter_policy(tmp_path):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"vcf\n")
    args = argparse.Namespace(
        modules=["formatter", "single_cell"],
        tools=["ldsc_celltype"],
        variant_id_type="unique",
        vcf=str(vcf),
        bcftools=sys.executable,
        output_directory=str(tmp_path / "results"),
        dataset_id="study",
        _step_num=1,
    )

    with pytest.raises(ValueError, match="HapMap3.*rsid"):
        run_formatter_runner(args, {})


def test_ldsc_celltype_pipeline_runner_injects_formatter_output_without_magma(
    tmp_path, monkeypatch,
):
    formatter_sumstats = tmp_path / "study_ldsc.tsv"
    formatter_sumstats.write_text("SNP\tZ\nrs1\t2.5\n", encoding="utf-8")
    args = argparse.Namespace(
        tools=["ldsc_celltype"],
        output_directory=str(tmp_path / "results"),
        dataset_id="study",
        _step_num=2,
    )
    context = {
        "formatter": {"ldsc": {"ldsc_file": str(formatter_sumstats)}},
    }
    observed = {}

    def fake_run(runtime_args, runtime_context):
        observed["source"] = runtime_args.ldsc_celltype_sumstats_source
        observed["sumstats"] = runtime_args.ldsc_celltype_sumstats_file
        observed["context"] = runtime_context
        return "completed"

    monkeypatch.setattr(
        "postgwas.modules.single_cell.service.run_single_cell_direct",
        fake_run,
    )

    result = run_single_cell_runner(args, context)

    assert result == "completed"
    assert observed["source"] == "formatter"
    assert observed["sumstats"] == str(formatter_sumstats)
    assert observed["context"] is context
    assert args.output_directory == str(tmp_path / "results")


@pytest.mark.skipif(
    os.environ.get("POSTGWAS_RUN_SCDRS_INTEGRATION") != "1",
    reason="set POSTGWAS_RUN_SCDRS_INTEGRATION=1 to run the pinned native CLI",
)
def test_scdrs_pinned_native_cli_smoke(tmp_path):
    anndata = pytest.importorskip("anndata")
    executable = shutil.which("scdrs")
    if executable is None:
        pytest.skip("pinned scDRS executable is not installed")
    random = np.random.default_rng(20220809)
    cells = ["cell_%03d" % index for index in range(80)]
    genes = ["GENE%03d" % index for index in range(400)]
    matrix = random.poisson(2, size=(len(cells), len(genes))).astype(float)
    adata = anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
        dtype=matrix.dtype,
    )
    h5ad = tmp_path / "native_atlas.h5ad"
    adata.write_h5ad(h5ad)
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_h5ad_file=str(h5ad),
        scdrs_gene_set_file=str(_gene_sets(tmp_path / "native.gs")),
        scdrs=str(executable),
        scdrs_control_gene_sets=5,
        dataset_id="native",
        output_directory=str(tmp_path / "results"),
    )

    result = run_single_cell_direct(args)

    assert result.artifacts["scdrs_score_1"].path.is_file()
    assert result.artifacts["scdrs_full_score_1"].path.is_file()


@pytest.mark.skipif(
    os.environ.get("POSTGWAS_RUN_SCDRS_INTEGRATION") != "1",
    reason="set POSTGWAS_RUN_SCDRS_INTEGRATION=1 to run the pinned native CLI",
)
def test_scdrs_pinned_native_magma_gene_set_smoke(tmp_path):
    anndata = pytest.importorskip("anndata")
    executable = shutil.which("scdrs")
    if executable is None:
        pytest.skip("pinned scDRS executable is not installed")
    random = np.random.default_rng(20220810)
    cells = ["cell_%03d" % index for index in range(80)]
    genes = ["GENE%03d" % index for index in range(5000)]
    matrix = random.poisson(2, size=(len(cells), len(genes))).astype(float)
    adata = anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
        dtype=matrix.dtype,
    )
    h5ad = tmp_path / "native_magma_atlas.h5ad"
    adata.write_h5ad(h5ad)
    args = argparse.Namespace(
        tools=["scdrs"],
        scdrs_gene_set_source="magma",
        scdrs_h5ad_file=str(h5ad),
        scdrs_magma_gene_results_file=str(
            _magma_gene_statistics(tmp_path / "native.genes.out")
        ),
        scdrs_gene_id_map=str(_identifier_map(tmp_path / "native_map.tsv")),
        scdrs=str(executable),
        scdrs_control_gene_sets=5,
        dataset_id="native_magma",
        output_directory=str(tmp_path / "results"),
    )

    result = run_single_cell_direct(args)

    assert result.artifacts["scdrs_generated_gene_set"].path.is_file()
    assert result.artifacts["scdrs_gene_mapping_report"].path.is_file()
    assert result.artifacts["scdrs_score_1"].path.is_file()
