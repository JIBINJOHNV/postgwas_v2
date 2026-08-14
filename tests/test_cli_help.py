"""Regression tests for the complete public command-line help surface."""

from __future__ import annotations

import argparse
from collections.abc import Collection
from io import StringIO
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

import pytest
from rich.console import Console
from rich.text import Text

from postgwas.config import load_configuration
from postgwas.modules.enrichment.cli import get_geneset_parser
from postgwas.modules.caldera.cli import build_parser as build_caldera_parser
from postgwas.modules.filtering.cli import build_parser as build_filtering_parser
from postgwas.modules.flames.cli import build_parser as build_flames_parser
from postgwas.modules.formatting.cli import build_parser as build_formatting_parser
from postgwas.modules.gcta_cojo.cli import build_parser as build_gcta_cojo_parser
from postgwas.modules.gcta_gene.cli import (
    build_parser as build_gcta_gene_parser,
    get_gcta_gene_parser,
)
from postgwas.modules.harmonisation.cli import get_harmonisation_parser
from postgwas.modules.harmonisation.concordance.cli import get_validation_parser
from postgwas.modules.imputation.cli import build_parser as build_imputation_parser
from postgwas.modules.ld_annotation.cli import (
    build_parser as build_ld_annotation_parser,
)
from postgwas.modules.ld_clumping.cli import build_parser as build_ld_clumping_parser
from postgwas.modules.ldsc.cli import build_parser as build_ldsc_parser
from postgwas.modules.magma.cli import build_parser as build_magma_parser
from postgwas.modules.magmacovar.cli import build_parser as build_magmacovar_parser
from postgwas.modules.manhattan.cli import build_parser as build_manhattan_parser
from postgwas.modules.mixer.cli import (
    build_parser as build_mixer_parser,
    get_mixer_parser,
)
from postgwas.modules.pops.cli import build_parser as build_pops_parser
from postgwas.modules.qc_summary.cli import build_parser as build_qc_parser
from postgwas.modules.single_cell.cli import build_parser as build_single_cell_parser
from postgwas.modules.fine_mapping.cli import build_parser as build_finemap_parser
from postgwas.modules.kpops import service as kpops_service
from postgwas.modules.kpops.cli import (
    build_parser as build_kpops_parser,
    get_kpops_direct_examples,
)
from postgwas.cli.common import (
    get_annot_ldblock_parser,
    get_assoc_plot_parser,
    get_common_imputation_parser,
    get_imputation_population_parser,
    get_ld_annotation_genome_build_parser,
    get_ld_clump_parser,
    get_ld_clumping_population_parser,
    get_pipeline_genome_build_parser,
    sumstat_summary_arg_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    format_cli_choices,
    format_cli_default,
    help_with_choices,
    help_with_conditional_requirement,
    help_with_default,
    mark_cli_required_help,
    style_cli_defaults,
    style_cli_requirement,
)
from postgwas.pipeline.registry import REGISTRY
from postgwas.pipeline import cli as pipeline_cli
from postgwas.core.errors import ConfigurationError


PUBLIC_HELP_COMMANDS = {
    "global": ("--help",),
    "annot_ldblock": ("annot_ldblock", "--help"),
    "config": ("config", "--help"),
    "config export": ("config", "export", "--help"),
    "config show": ("config", "show", "--help"),
    "config validate": ("config", "validate", "--help"),
    "caldera": ("caldera", "--help"),
    "finemap": ("finemap", "--help"),
    "flames": ("flames", "--help"),
    "formatter": ("formatter", "--help"),
    "gcta_gene": ("gcta_gene", "--help"),
    "gcta_cojo": ("gcta_cojo", "--help"),
    "harmonisation": ("harmonisation", "--help"),
    "heritability": ("heritability", "--help"),
    "imputation": ("imputation", "--help"),
    "ld_clump": ("ld_clump", "--help"),
    "magma": ("magma", "--help"),
    "magmacovar": ("magmacovar", "--help"),
    "manhattan": ("manhattan", "--help"),
    "mixer": ("mixer", "--help"),
    "kpops": ("kpops", "--help"),
    "pathway_enrichment": ("pathway_enrichment", "--help"),
    "pipeline": ("pipeline", "--help"),
    "pipeline module": ("pipeline", "--modules", "mixer", "--help"),
    "pops": ("pops", "--help"),
    "qc": ("qc", "--help"),
    "resources": ("resources", "--help"),
    "single_cell": ("single_cell", "--help"),
    "sumstat_filter": ("sumstat_filter", "--help"),
    "validate": ("--validate", "--help"),
}

PIPELINE_TARGETS = tuple(
    name for name in REGISTRY.names() if REGISTRY.get(name).pipeline_enabled
)

DIRECT_ANALYSIS_PARSER_BUILDERS = {
    "annot_ldblock": build_ld_annotation_parser,
    "caldera": build_caldera_parser,
    "finemap": build_finemap_parser,
    "flames": build_flames_parser,
    "formatter": build_formatting_parser,
    "gcta_cojo": build_gcta_cojo_parser,
    "gcta_gene": build_gcta_gene_parser,
    "harmonisation": lambda: get_harmonisation_parser(add_help=True),
    "heritability": build_ldsc_parser,
    "imputation": build_imputation_parser,
    "kpops": build_kpops_parser,
    "ld_clump": build_ld_clumping_parser,
    "magma": build_magma_parser,
    "magmacovar": build_magmacovar_parser,
    "manhattan": build_manhattan_parser,
    "mixer": build_mixer_parser,
    "pathway_enrichment": lambda: get_geneset_parser(add_help=True),
    "pops": build_pops_parser,
    "qc": build_qc_parser,
    "single_cell": build_single_cell_parser,
    "sumstat_filter": build_filtering_parser,
    "validate": get_validation_parser,
}


def _help(arguments: tuple[str, ...], width: int) -> str:
    environment = os.environ.copy()
    environment.update({"COLUMNS": str(width), "PYTHONDONTWRITEBYTECODE": "1"})
    completed = subprocess.run(
        [sys.executable, "-m", "postgwas", *arguments],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout
    return completed.stdout


def _assert_only_standard_default_labels(help_text: str, context: str) -> None:
    """Require the exact unqualified, case-sensitive ``Default:`` label."""
    for match in re.finditer(r"\bdefault\s*:", help_text, flags=re.IGNORECASE):
        assert match.group(0) == "Default:", context
    qualified = re.search(
        r"\b(?:[A-Za-z][A-Za-z-]*[ \t]+)+default\s*:",
        help_text,
        flags=re.IGNORECASE,
    )
    assert qualified is None, (context, qualified.group(0) if qualified else None)


def _visible_finite_choice_actions(parser: argparse.ArgumentParser):
    actions = []
    for action in parser._actions:
        if not isinstance(action.choices, Collection) or not action.choices:
            continue
        if action.help is argparse.SUPPRESS:
            continue
        assert action.help is not None, action.dest
        actions.append(action)
    return actions


def _assert_all_finite_choices_are_visible(parser, context: str) -> None:
    normalized_help = " ".join(parser.format_help().split())
    actions = _visible_finite_choice_actions(parser)
    assert normalized_help.count("Available options:") == len(actions), context
    for action in actions:
        rendered = Text.from_markup(format_cli_choices(action.choices)).plain
        assert rendered in normalized_help, (context, action.dest)


def test_public_help_matrix_covers_every_registered_command():
    covered_commands = {
        arguments[0]
        for arguments in PUBLIC_HELP_COMMANDS.values()
        if not arguments[0].startswith("-")
    }
    assert covered_commands == set(REGISTRY.commands())


def test_direct_choice_audit_covers_every_public_analysis_command():
    expected = (
        set(REGISTRY.commands()) - {"config", "pipeline", "resources"}
    ) | {"validate"}
    assert set(DIRECT_ANALYSIS_PARSER_BUILDERS) == expected

    for command, builder in DIRECT_ANALYSIS_PARSER_BUILDERS.items():
        _assert_all_finite_choices_are_visible(builder(), command)


@pytest.mark.parametrize("width", (80, 120))
def test_every_public_help_page_is_readable_and_has_examples(width):
    for command, arguments in PUBLIC_HELP_COMMANDS.items():
        help_text = _help(arguments, width)
        assert "postgwas" in help_text, command
        assert "example" in help_text.lower(), command
        assert "usage: postgwas-" not in help_text.lower(), command
        assert "[cyan]" not in help_text, command
        assert "[/cyan]" not in help_text, command
        assert "[/bold green]" not in help_text, command
        assert "REQUIRED." not in help_text, command
        assert "scientific" not in help_text.lower(), command
        _assert_only_standard_default_labels(help_text, command)
        if command not in {
            "config",
            "config export",
            "config show",
            "config validate",
            "resources",
        }:
            assert "--show-screen" in help_text, command
            assert "--hide-screen" in help_text, command
            assert "--resume" in help_text, command
            assert "--no-resume" in help_text, command
            assert "--overwrite" in help_text, command
            assert "Default: true" in help_text, command
            assert "Default: false" in help_text, command


@pytest.mark.parametrize("width", (80, 120))
@pytest.mark.parametrize(
    ("arguments", "build_marker", "population_marker"),
    (
        (("gcta_gene", "--help"), "Required:", "Required:"),
        (
            ("pipeline", "--modules", "gcta_gene", "--help"),
            "Shared genome build",
            "Population represented",
        ),
    ),
)
def test_cli_help_uses_shared_description_column_for_same_and_wrapped_rows(
    width, arguments, build_marker, population_marker,
):
    help_lines = _help(arguments, width).splitlines()
    build_line = next(
        line for line in help_lines if "--genome-build BUILD" in line
    )
    build_index = help_lines.index(build_line)
    description_column = build_line.index(build_marker)

    assert build_line.lstrip().startswith("--genome-build BUILD")
    assert build_marker in build_line
    assert (
        len(help_lines[build_index + 1])
        - len(help_lines[build_index + 1].lstrip())
        == description_column
    )

    population_index = next(
        index
        for index, line in enumerate(help_lines)
        if "--gcta-reference-population CODE" in line
    )
    population_description = next(
        line
        for line in help_lines[population_index:population_index + 2]
        if population_marker in line
    )
    if width == 80:
        assert population_marker not in help_lines[population_index]
        assert population_description == help_lines[population_index + 1]
    else:
        assert population_marker in help_lines[population_index]
    assert description_column == population_description.index(
        population_marker
    )


def test_shared_runtime_flags_leave_yaml_authoritative():
    parser = get_compute_parser()

    assert not hasattr(parser.parse_args([]), "show_screen")
    assert not hasattr(parser.parse_args([]), "resume")
    assert not hasattr(parser.parse_args([]), "overwrite")
    assert parser.parse_args(["--show-screen"]).show_screen is True
    assert parser.parse_args(["--hide-screen"]).show_screen is False
    assert parser.parse_args(["--resume"]).resume is True
    assert parser.parse_args(["--no-resume"]).resume is False
    assert parser.parse_args(["--overwrite"]).overwrite is True
    with pytest.raises(SystemExit):
        parser.parse_args(["--show-screen", "--hide-screen"])

    actions = {
        action.dest: action
        for action in parser._actions
        if action.dest in {"resume", "overwrite"}
    }
    assert "[bold green]Default:[/bold green] [green]true[/green]" in (
        actions["resume"].help or ""
    )
    assert "[bold green]Default:[/bold green] [green]false[/green]" in (
        actions["overwrite"].help or ""
    )


def test_every_pipeline_target_builds_context_specific_help():
    for module in PIPELINE_TARGETS:
        help_text = _help(("pipeline", "--modules", module, "--help"), 120)
        assert "Usage: postgwas pipeline [options] --modules %s" % module in help_text
        assert "Workflow examples" in help_text
        assert "--show-screen" in help_text
        assert "--hide-screen" in help_text
        assert "--resume" in help_text
        assert "--no-resume" in help_text
        assert "--overwrite" in help_text
        assert "Default: true" in help_text
        assert "Default: false" in help_text
        assert "scientific" not in help_text.lower(), module
        _assert_only_standard_default_labels(help_text, "pipeline:%s" % module)


def test_every_pipeline_target_displays_every_finite_choice(monkeypatch):
    captured = {}

    def capture_help(modules, parser, *, error=False):
        captured["parser"] = parser

    monkeypatch.setattr(pipeline_cli, "print_full_pipeline_help", capture_help)
    for module in PIPELINE_TARGETS:
        captured.clear()
        monkeypatch.setattr(
            sys,
            "argv",
            ["postgwas pipeline", "--modules", module, "--help"],
        )
        with pytest.raises(SystemExit) as stopped:
            pipeline_cli.main()
        assert stopped.value.code == 0, module
        parser = captured["parser"]
        modules_action = next(
            action for action in parser._actions if action.dest == "modules"
        )
        assert tuple(modules_action.choices or ()) == REGISTRY.names()
        _assert_all_finite_choices_are_visible(parser, "pipeline:%s" % module)
        _assert_only_standard_default_labels(
            parser.format_help(), "pipeline:%s" % module,
        )


def test_generic_pipeline_help_lists_modules_and_single_cell_tools():
    help_text = " ".join(_help(("pipeline", "--help"), 160).split())
    expected_modules = Text.from_markup(
        format_cli_choices(REGISTRY.names())
    ).plain
    expected_tools = Text.from_markup(
        format_cli_choices(pipeline_cli.SINGLE_CELL_TOOLS)
    ).plain

    assert expected_modules in help_text
    assert expected_tools in help_text


def test_readme_prose_avoids_unnecessary_scientific_wording():
    repository = Path(__file__).parents[1]
    for readme in repository.rglob("*.md"):
        relative = readme.relative_to(repository)
        if "readme" not in readme.name.lower() or any(
            part.startswith(".") or part == "_to_delete"
            for part in relative.parts[:-1]
        ):
            continue
        visible_text = re.sub(
            r"\]\([^)]*\)", "]", readme.read_text(encoding="utf-8"),
        )
        assert "scientific" not in visible_text.lower(), readme


def test_pipeline_magma_help_hides_only_formatter_created_inputs():
    help_text = _help(("pipeline", "--modules", "magma", "--help"), 120)
    normalized_help = " ".join(help_text.split())

    assert "--resume, --no-resume" in help_text
    assert "--overwrite" in help_text
    assert "--resolve-variants-to-reference" in help_text
    assert "Default: false" in normalized_help
    assert "Default: true" in help_text
    assert "--snp-location-file" not in help_text
    assert "--p-value-file" not in help_text
    assert "--variant-id-type" not in help_text
    for external_resource in (
        "--magma-ld-reference",
        "--gene-location-file",
        "--gene-set-file",
    ):
        assert external_resource in help_text
    assert "--magma-ld-reference PREFIX Required:" in normalized_help
    assert "--gene-location-file PATH Required:" in normalized_help
    assert "--gene-set-file PATH Required:" not in normalized_help
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("magma").required_options
    }
    assert required["magma_ld_reference"] == (
        "modules.magma.input.ld_reference_prefix"
    )
    assert required["gene_location_file"] == (
        "modules.magma.input.gene_location_file"
    )


def test_pipeline_magma_required_references_may_come_from_yaml(
    tmp_path, monkeypatch,
):
    vcf = tmp_path / "study.vcf.gz"
    vcf.write_bytes(b"placeholder")
    run_config = tmp_path / "magma.yaml"
    run_config.write_text(
        "modules:\n"
        "  magma:\n"
        "    input:\n"
        "      ld_reference_prefix: /configured/reference\n"
        "      gene_location_file: /configured/genes.loc\n",
        encoding="utf-8",
    )
    executed = []
    monkeypatch.setattr(
        pipeline_cli,
        "execute_pipeline",
        lambda args, plan, configuration: executed.append(
            (args, plan, configuration)
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "postgwas pipeline",
            "--modules",
            "magma",
            "--vcf",
            str(vcf),
            "--dataset-id",
            "STUDY",
            "--output-directory",
            str(tmp_path / "results"),
            "--run-config",
            str(run_config),
        ],
    )

    assert pipeline_cli.main() is None
    assert len(executed) == 1


def test_pipeline_gcta_help_shows_analysis_options_and_specific_examples():
    help_text = _help(("pipeline", "--modules", "gcta_gene", "--help"), 120)
    normalized_help = " ".join(help_text.split())

    assert "--gcta-input-file" not in help_text
    for option in (
        "--method",
        "--gcta-reference-prefix",
        "--gene-list",
        "--fastbat-set-list",
        "--gmt",
        "--genome-build",
        "--gcta-reference-population",
    ):
        assert option in help_text
    for description in (
        "fastbat_gene: Test SNPs assigned to each gene and its configured window",
        "fastbat_segment: Test consecutive fixed-size genomic segments",
        "fastbat_set: Test custom SNP sets supplied through --fastbat-set-list",
        "mbat_combo: Test each gene by combining signed mBAT and unsigned fastBAT",
    ):
        assert description in normalized_help
        assert any(
            line.lstrip().startswith(description.split(":", 1)[0] + ":")
            for line in help_text.splitlines()
        )
    assert "Run gene-based mBAT-combo:" in help_text
    assert "--method mbat_combo" in help_text
    assert "Convert a GMT pathway file and run custom-set fastBAT:" in help_text
    assert "--method fastbat_set" in help_text
    assert "--gmt pathways.gmt" in help_text


def test_pipeline_pops_help_shows_required_external_resources():
    direct_help = build_pops_parser().format_help()
    help_text = _help(("pipeline", "--modules", "pops", "--help"), 120)

    assert "PoPS does not use MAGMACOVAR .gsa.out results" in direct_help
    assert "technical covariates read from MAGMA .genes.raw" in direct_help
    assert "--gene-universe-policy POLICY" in direct_help
    assert "--gene-universe-policy intersect" in direct_help
    assert "Run MAGMA gene analysis followed by PoPS prioritisation:" in help_text
    assert "Run MAGMA and explicitly intersect incompatible PoPS target genes:" in help_text
    assert "--gene-universe-policy intersect" in help_text
    for option in (
        "--magma-ld-reference reference/1000G_EUR",
        "--gene-location-file reference/FUMA/ENSGv102.coding.genes.txt",
        "--feature-matrix-prefix reference/FLAMES/pops_features_full_FUMA_compatible/features_munged/pops_features",
        "--feature-matrix-chunks 116",
        "--pops-gene-location-file reference/FLAMES/pops_features_full_FUMA_compatible/gene_annots.txt",
        "--control-features-file reference/FLAMES/pops_features_full_FUMA_compatible/control.features",
        "--pops-genome-build GRCh37",
    ):
        assert option in help_text
    assert "--run-config analysis_pipeline.yaml" not in help_text


def test_finemap_help_has_runnable_direct_and_pipeline_examples():
    direct_help = _help(("finemap", "--help"), 120)
    pipeline_help = _help(("pipeline", "--modules", "finemap", "--help"), 120)

    assert "--window-kb INT" in direct_help
    assert "--window-kb INT" in pipeline_help
    assert "--ld-window-kb KB" in pipeline_help

    for option in (
        "--susie-input-file formatted/STUDY_susie.tsv.gz",
        "--locus-file loci.tsv",
        "--finemap-ld-reference reference/1000G_EUR",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in direct_help
    for option in (
        "--vcf study_GRCh37.vcf.gz",
        "--genome-build GRCh37",
        "--ld-region-dir reference/ld_blocks",
        "--ld-block-populations EUR",
        "--ld-folder reference/pairwise_ld",
        "--population EUR",
        "--variant-id-type unique",
        "--finemap-method susie",
        "--finemap-ld-reference reference/1000G_EUR",
        "--plink plink",
        "--bcftools bcftools",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help


def test_kpops_help_has_runnable_direct_and_dependency_complete_pipeline_examples():
    direct_help = build_kpops_parser().format_help()
    pipeline_help = _help(("pipeline", "--modules", "kpops", "--help"), 120)
    normalized_pipeline_help = " ".join(pipeline_help.split())

    for option in (
        (
            "--magma-association-prefix "
            "magma/02_intermediates/positional/native_outputs/"
            "STUDY_magma_35up_10down"
        ),
        (
            "--kpops-gene-annotation-file /path/to/postgwas-resources/"
            "pops/GRCh37_gene_annot_jun10.txt"
        ),
        (
            "--kernel-matrix-prefix /path/to/postgwas-resources/kpops/"
            "kernels/GRCh37/pops_features_standardized_linear"
        ),
        "--kpops-genome-build GRCh37",
    ):
        assert option in direct_help
    assert "--gene-universe-policy POLICY" in direct_help
    assert "Default: intersect" in direct_help
    for option in (
        "--vcf study_GRCh37.vcf.gz",
        (
            "--magma-ld-reference /path/to/postgwas-resources/magma/"
            "functional_mapping/base/ld_reference/g1000_eur/g1000_eur"
        ),
        (
            "--gene-location-file /path/to/postgwas-resources/kpops/"
            "software/8acd49ed8c96565b17c2997420f514f24b65e096/data/"
            "Ensembl.hg19.gene.loc"
        ),
        (
            "--kpops-gene-annotation-file /path/to/postgwas-resources/"
            "pops/GRCh37_gene_annot_jun10.txt"
        ),
        (
            "--kernel-matrix-prefix /path/to/postgwas-resources/kpops/"
            "kernels/GRCh37/pops_features_standardized_linear"
        ),
        "--kpops-genome-build GRCh37",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--gene-universe-policy POLICY" in pipeline_help
    for option in (
        "--kpops-gene-annotation-file PATH",
        "--kernel-matrix-prefix PREFIX",
        "--kpops-genome-build BUILD",
    ):
        assert "%s Required:" % option in normalized_pipeline_help
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("kpops").required_options
    }
    assert required["kpops_gene_annotation_file"] == (
        "modules.kpops.gene_annotation_file"
    )
    assert required["kernel_matrix_prefix"] == (
        "modules.kpops.kernel_matrix_prefix"
    )
    assert required["kpops_genome_build"] == "modules.kpops.genome_build"
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help


def test_kpops_upstream_apob_example_is_parseable_and_preserves_model_settings():
    examples = get_kpops_direct_examples()
    title, command, arguments = examples[1]
    assert title == "Reproduce the upstream v1.0.0 ApoB command settings:"
    assert command == "postgwas kpops"

    parsed = build_kpops_parser().parse_args(
        [token for argument in arguments for token in shlex.split(argument)]
    )

    assert parsed.magma_association_prefix == "/path/to/ApoB"
    assert parsed.kpops_gene_annotation_file.endswith(
        "/pops/GRCh37_gene_annot_jun10.txt"
    )
    assert parsed.kernel_matrix_prefix == "/path/to/ApoB/kernel_linear"
    assert parsed.kpops_genome_build == "GRCh37"
    assert parsed.gene_universe_policy == "strict"
    assert parsed.use_magma_covariates is True
    assert parsed.training_chromosomes == ["loco"]
    assert parsed.kpops_device == "cuda"
    assert parsed.top_contributor_gene_count == 5
    assert parsed.anchor_genes == [
        "PCSK9", "ANGPTL3", "APOB", "ABCG5", "ALB", "NPC1L1", "GIGYF1",
        "JAK2", "A1CF", "PDE3B", "APOC3", "CETP", "ASGR1", "LDLR",
        "ZNF234", "ZNF229", "NECTIN2", "RRBP1",
    ]
    assert parsed.anchor_gene_type == "NAME"
    assert parsed.seed == "42"
    assert parsed.kpops_verbose is True
    assert parsed.dataset_id == "ApoB"
    assert parsed.output_directory == "results/ApoB"

    configuration = kpops_service._resolved_configuration(parsed)
    upstream = kpops_service._arguments(
        configuration,
        Path("results/ApoB/ApoB_kpops"),
        Path("k-pops.py"),
    )
    assert upstream[upstream.index("--random_seed") + 1] == "42"
    assert upstream[upstream.index("--device") + 1] == "cuda"
    assert upstream[upstream.index("--top_n_contributor_gene") + 1] == "5"
    assert upstream[upstream.index("--training_chromosomes") + 1] == "loco"
    assert "--use_magma_covariates" in upstream
    assert "--project_out_covariates_remove_hla" in upstream
    assert "--training_remove_hla" in upstream
    assert "--testing_remove_hla" in upstream
    assert "--no_save_attr_files" in upstream
    assert "--verbose" in upstream
    assert upstream[upstream.index("--anchor_genes_type") + 1] == "NAME"
    assert upstream[upstream.index("--anchor_genes") + 1] == (
        "PCSK9,ANGPTL3,APOB,ABCG5,ALB,NPC1L1,GIGYF1,JAK2,A1CF,PDE3B,"
        "APOC3,CETP,ASGR1,LDLR,ZNF234,ZNF229,NECTIN2,RRBP1"
    )


def test_kpops_help_explains_every_public_argument_and_model_effect():
    parser = build_kpops_parser()
    help_text = " ".join(parser.format_help().split())

    for action in parser._actions:
        if action.dest != "help":
            assert action.help not in (None, argparse.SUPPRESS), action.dest
    for explanation in (
        "Both files must come from the same gene-association run",
        "float32 matrix must be square, finite, symmetric",
        "Installation normally provides it",
        "does not infer or lift coordinates",
        "leave-one-chromosome-out model",
        "predict the remaining chromosomes",
        "cuda for a supported NVIDIA GPU",
        "changes explanation columns, not the K-POPS score",
        "calculate each prediction's Anchor_Score",
        "Ambiguous annotation symbols are rejected",
        "derived from the .genes.raw file out of target Z statistics",
        "Storage grows quadratically with gene count",
        "captured command log",
        "command-line values override matching YAML keys",
    ):
        assert explanation in help_text
    for default_value in (
        "loco",
        "cpu",
        "5",
        "none",
        "ENSGID",
        "true",
        "false",
        "k-pops.py",
    ):
        assert "Default: %s" % default_value in help_text


def test_caldera_help_has_runnable_direct_and_complete_dependency_pipeline_examples():
    direct_help = build_caldera_parser().format_help()
    pipeline_help = _help(("pipeline", "--modules", "caldera", "--help"), 120)

    for option in (
        "--pops-file pops/STUDY_pops.preds",
        "--credible-set-file finemap/STUDY_credible_sets_GRCh37.tsv",
        "--caldera-genome-build GRCh37",
    ):
        assert option in direct_help
    for option in (
        "--vcf study_GRCh37.vcf.gz",
        "--genome-build GRCh37",
        "--ld-region-dir reference/ld_blocks",
        "--ld-block-populations EUR",
        "--ld-folder reference/ld",
        "--population EUR",
        "--magma-ld-reference reference/1000G_EUR",
        "--gene-location-file reference/NCBI37.3.gene.loc",
        "--feature-matrix-prefix reference/pops/features_munged/pops_features",
        "--feature-matrix-chunks 116",
        "--pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv",
        "--pops-genome-build GRCh37",
        "--finemap-method susie",
        "--finemap-ld-reference reference/1000G_EUR",
        "--caldera-genome-build GRCh37",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help


def test_flames_help_has_runnable_direct_and_dependency_complete_pipeline_examples():
    direct_help = build_flames_parser().format_help()
    pipeline_help = _help(("pipeline", "--modules", "flames", "--help"), 120)

    assert direct_help.count("--tabix PATH") == 1
    assert pipeline_help.count("--tabix PATH") == 1

    for option in (
        "--credible-sets-directory finemap/downstream_inputs/flames",
        "--magma-gene-results-file magma/STUDY.genes.out",
        "--magma-covariate-results-file magma/STUDY.gsa.out",
        "--pops-scores-file pops/STUDY.preds",
        "--flames-annotation-directory reference/FLAMES/Annotation_data",
        "--flames-genome-build GRCh37",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in direct_help
    for option in (
        "--modules flames",
        "--vcf study_GRCh37.vcf.gz",
        "--genome-build GRCh37",
        "--ld-region-dir reference/ld_blocks",
        "--ld-block-populations EUR",
        "--ld-folder reference/pairwise_ld",
        "--population EUR",
        "--magma-ld-reference reference/1000G_EUR",
        "--gene-location-file reference/FUMA/ENSGv102.coding.genes.txt",
        "--covariates reference/GTEx/gtex_v8_ts_avg_log2TPM.txt",
        "--feature-matrix-prefix reference/FLAMES/pops_features_full_FUMA_compatible/features_munged/pops_features",
        "--feature-matrix-chunks 116",
        "--pops-gene-location-file reference/FLAMES/pops_features_full_FUMA_compatible/gene_annots.txt",
        "--control-features-file reference/FLAMES/pops_features_full_FUMA_compatible/control.features",
        "--pops-genome-build GRCh37",
        "--finemap-method susie",
        "--finemap-ld-reference reference/1000G_EUR",
        "--flames-annotation-directory reference/FLAMES/Annotation_data",
        "--flames-genome-build GRCh37",
        "--plink plink",
        "--bcftools bcftools",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help


def test_magmacovar_help_distinguishes_marginal_and_tissue_specific_models():
    configured = load_configuration().modules.magmacovar
    parser = build_magmacovar_parser()
    help_text = parser.format_help()
    normalized_help = " ".join(help_text.split())

    model_action = next(
        action for action in parser._actions if action.dest == "covariate_model"
    )
    direction_action = next(
        action for action in parser._actions if action.dest == "covariate_direction"
    )
    assert "Accepted modifiers: %s" % ", ".join(configured.model_modifiers) in (
        model_action.help or ""
    )
    assert tuple(direction_action.choices or ()) == (
        "two-sided",
        "greater",
        "smaller",
    )
    assert "--covariate-direction {two-sided,greater,smaller}" in help_text

    for placeholder, description in configured.model_placeholders.items():
        assert "<%s>: %s" % (placeholder, description) in normalized_help
    for specification in configured.model_modifiers.values():
        for syntax in specification.syntax:
            assert syntax in normalized_help
        assert specification.description in normalized_help
        assert specification.published_example in normalized_help
        assert specification.example in normalized_help
    assert "How to choose a MAGMA model" in help_text
    assert "Question: Does each property associate" in normalized_help
    assert "Use: omit --covariate-model" in normalized_help
    assert "Direction: --covariate-direction two-sided" in normalized_help
    assert "Question: Is a tissue positively associated" in normalized_help
    assert "Use: --covariate-model condition-hide=Average" in normalized_help
    assert "Direction: --covariate-direction greater" in normalized_help
    assert "Advanced MAGMA model modifiers" in help_text
    assert "Example property names are illustrative" in normalized_help
    assert normalized_help.count("Published GWAS example:") == len(
        configured.model_modifiers
    )
    assert (
        "Interaction and gene-set-only modifiers are unavailable" in normalized_help
    )
    assert "Omit this option to test every property separately" in normalized_help
    assert "Two-sided is MAGMA's default for any gene property" in normalized_help
    assert "original FLAMES README" in normalized_help
    assert (
        "Create tissue-relevance input as documented by original FLAMES"
        in normalized_help
    )
    assert "gtex_v8_ts_avg_log2TPM.txt" in normalized_help
    assert "separate FUMA-compatible tissue-specific model" in normalized_help
    assert "gtex_v8_ts_general_avg_log2TPM.txt" in normalized_help
    assert "--covariate-model condition-hide=Average" in normalized_help
    assert "--covariate-direction greater" in normalized_help
    assert "54-specific-tissue" in normalized_help
    assert "does not match the literal command printed in its README" in normalized_help


def test_direct_magma_help_retains_direct_input_files():
    help_text = _help(("magma", "--help"), 120)
    normalized_help = " ".join(help_text.split())
    assert "--snp-location-file" in help_text
    assert "--p-value-file" in help_text
    assert "--resolve-variants-to-reference" in help_text
    assert "Default: false" in normalized_help
    for required_option in (
        "--snp-location-file PATH",
        "--p-value-file PATH",
        "--magma-ld-reference PREFIX",
        "--gene-location-file PATH",
        "--dataset-id NAME",
        "--output-directory PATH",
    ):
        assert "%s Required:" % required_option in normalized_help
    assert "--gene-set-file PATH Required:" not in normalized_help


def test_long_option_help_wraps_and_aligns_at_narrow_terminal_width():
    help_text = _help(("pipeline", "--modules", "magma", "--help"), 80)
    lines = help_text.splitlines()
    memory_index = next(
        index for index, line in enumerate(lines) if "--memory-gb GB" in line
    )
    memory_line = lines[memory_index]
    description_column = memory_line.index("Maximum total memory")

    assert len(memory_line) <= 80
    assert lines[memory_index + 1].startswith(" " * description_column)
    assert lines[memory_index + 1].lstrip().startswith("If omitted, PostGWAS uses")


def test_aligned_help_preserves_deliberate_option_rows(monkeypatch):
    monkeypatch.setenv("COLUMNS", "100")
    parser = argparse.ArgumentParser(formatter_class=AlignedRichHelpFormatter)
    parser.add_argument(
        "--method",
        help="Choose a method:\nalpha: First analysis.\nbeta: Second analysis.",
    )

    lines = [line.strip() for line in parser.format_help().splitlines()]

    assert any(line.endswith("Choose a method:") for line in lines)
    assert "alpha: First analysis." in lines
    assert "beta: Second analysis." in lines


@pytest.mark.parametrize(
    ("builder", "destinations"),
    (
        (build_magma_parser, {"snp_location_file", "p_value_file"}),
        (build_magmacovar_parser, {"magma_gene_results_file"}),
        (build_pops_parser, {"magma_association_prefix", "verbose"}),
        (lambda: get_geneset_parser(add_help=True), {"gene_input_file"}),
    ),
)
def test_user_friendly_flags_preserve_runner_destinations(builder, destinations):
    parser = builder()
    parser_destinations = {action.dest for action in parser._actions}
    assert destinations <= parser_destinations


def test_default_label_and_value_are_green_through_central_styles():
    markup = format_cli_default(50)
    assert markup == (
        "[bold green]Default:[/bold green] [green]50[/green]"
    )
    rendered = Text.from_markup(markup)
    assert rendered.plain == "Default: 50"
    assert [str(span.style) for span in rendered.spans] == ["bold green", "green"]
    with pytest.raises(TypeError):
        format_cli_default(False, label="Configured default")
    assert help_with_default("Window size", "[configured]") == (
        "Window size ([bold green]Default:[/bold green] "
        "[green]\\[configured][/green])."
    )
    assert style_cli_defaults("Threshold (default: 0.2).") == (
        "Threshold ([bold green]Default:[/bold green] [green]0.2[/green])."
    )
    assert style_cli_defaults("Packaged default: EUR") == (
        "[bold green]Default:[/bold green] [green]EUR[/green]"
    )
    assert style_cli_defaults("Configured YAML default: false") == (
        "[bold green]Default:[/bold green] [green]false[/green]"
    )
    legacy_styled = "[bold green]Default:[/bold green] [cyan]0.2[/cyan]"
    assert style_cli_defaults(legacy_styled) == (
        "[bold green]Default:[/bold green] [green]0.2[/green]"
    )
    assert style_cli_defaults("Configured default: [cyan]sss[/cyan].") == (
        "[bold green]Default:[/bold green] [green]sss[/green]."
    )
    already_styled = "[bold green]Default:[/bold green] [green]0.2[/green]"
    assert style_cli_defaults(already_styled) == already_styled
    qualified_styled = (
        "[bold green]Configured default:[/bold green] [green]0.2[/green]"
    )
    assert style_cli_defaults(qualified_styled) == already_styled


def test_public_parsers_do_not_define_or_parameterize_default_labels():
    repository = Path(__file__).parents[1]
    module_root = repository / "src" / "postgwas" / "modules"
    public_help_sources = {
        repository / "src" / "postgwas" / "cli" / "common.py",
        repository / "src" / "postgwas" / "cli" / "compute.py",
        repository / "src" / "postgwas" / "config" / "cli.py",
        repository / "src" / "postgwas" / "pipeline" / "cli.py",
        repository / "src" / "postgwas" / "resources" / "cli.py",
        module_root / "enrichment" / "arguments.py",
        module_root / "fine_mapping" / "arguments.py",
        *module_root.rglob("cli.py"),
    }
    forbidden = (
        re.compile(r"label\s*=\s*['\"][^'\"]*default", re.IGNORECASE),
        re.compile(r"\[bold green\][^\n]*default\s*:", re.IGNORECASE),
        re.compile(
            r"(?:configured(?: YAML)?|effective|packaged|upstream) "
            r"default\s*:",
            re.IGNORECASE,
        ),
    )

    for path in public_help_sources:
        source = path.read_text(encoding="utf-8")
        for pattern in forbidden:
            assert pattern.search(source) is None, (path, pattern.pattern)


def test_all_finite_choices_are_bright_yellow_without_changing_validation():
    markup = format_cli_choices(("cpu", "cuda", "mps"))
    assert markup == (
        "[bold bright_yellow]Available options:[/bold bright_yellow] "
        "[bright_yellow]cpu, cuda, mps[/bright_yellow]"
    )
    rendered = Text.from_markup(markup)
    assert rendered.plain == "Available options: cpu, cuda, mps"
    assert [str(span.style) for span in rendered.spans] == [
        "bold bright_yellow", "bright_yellow",
    ]

    parser = argparse.ArgumentParser(formatter_class=AlignedRichHelpFormatter)
    parser.add_argument(
        "--device",
        choices=("cpu", "cuda", "mps"),
        metavar="DEVICE",
        help="Compute device.",
    )
    parser.add_argument(
        "--tools",
        choices=("magma", "scdrs"),
        nargs="+",
        metavar="TOOL",
        help="One or more tools.",
    )
    parser.add_argument(
        "--minimum",
        choices=("5%", "10%"),
        metavar="PERCENT",
        help="Minimum percentage.",
    )
    parser.add_argument(
        "--backend",
        choices=("local", "container"),
        metavar="BACKEND",
        help=help_with_choices(
            "Execution backend.", ("local", "container"),
        ),
    )

    help_text = parser.format_help()
    assert "Available options: cpu, cuda, mps" in help_text
    assert "Available options: magma, scdrs" in help_text
    assert "Available options: 5%, 10%" in help_text
    normalized_help = " ".join(help_text.split())
    assert normalized_help.count("Available options: local, container") == 1
    assert parser.parse_args(["--device", "mps"]).device == "mps"
    with pytest.raises(SystemExit):
        parser.parse_args(["--device", "invalid"])


def test_formatter_multivalue_formats_are_prominent_and_complete():
    configuration = load_configuration()
    choices = tuple(configuration.modules.formatting.format_order)
    help_text = " ".join(_help(("formatter", "--help"), 160).split())
    expected = Text.from_markup(format_cli_choices(choices)).plain

    assert expected in help_text
    assert help_text.index(expected) < help_text.index("Format details:")
    for target in choices:
        assert target in help_text


def test_direct_and_pipeline_help_share_available_option_rendering():
    direct_help = _help(("kpops", "--help"), 120)
    pipeline_help = _help(("pipeline", "--modules", "kpops", "--help"), 120)

    for help_text in (direct_help, pipeline_help):
        normalized = " ".join(help_text.split())
        assert "Available options: strict, intersect" in normalized
        assert "Available options: GRCh37, GRCh38" in normalized
        assert "Available options: cpu, cuda, mps" in normalized
        assert "Available options: ENSGID, NAME" in normalized


def test_every_displayed_module_default_renders_entirely_in_green():
    builders = (
        get_compute_parser,
        lambda: get_gcta_gene_parser(direct_controls=True),
        get_harmonisation_parser,
        lambda: get_mixer_parser(direct_controls=True),
        build_magma_parser,
    )
    displayed_defaults = 0
    for builder in builders:
        for action in builder()._actions:
            help_text = action.help or ""
            if "default:" not in help_text.lower():
                continue
            displayed_defaults += 1
            rendered_help = style_cli_defaults(help_text)
            assert "[bold green]" in rendered_help
            assert "[green]" in rendered_help
            assert "[cyan]" not in rendered_help
    assert displayed_defaults >= 20


def test_legacy_module_parser_defaults_are_loaded_from_canonical_yaml():
    configuration = load_configuration()

    def action_default(parser, destination):
        return next(
            action.default
            for action in parser._actions
            if action.dest == destination
        )

    qc = configuration.modules.qc_summary
    qc_parser = sumstat_summary_arg_parser()
    assert action_default(qc_parser, "reference_af_column") == argparse.SUPPRESS
    assert action_default(qc_parser, "maximum_af_difference") == argparse.SUPPRESS
    assert qc.rules.maximum_af_difference == 0.2

    annotation = configuration.modules.ld_annotation
    assert action_default(
        get_ld_annotation_genome_build_parser(), "genome_build",
    ) == annotation.genome_build.value
    assert action_default(
        get_annot_ldblock_parser(), "ld_block_populations",
    ) == [population.value for population in annotation.populations]

    clumping = configuration.modules.ld_clumping
    clumping_parser = get_ld_clump_parser()
    for destination in (
        "clumping_methods", "lead_p", "candidate_p", "r2_clump", "r2_lead",
        "window_kb", "merge_dist", "missing_index_action", "remove_mhc",
    ):
        assert action_default(clumping_parser, destination) == argparse.SUPPRESS
    assert action_default(
        get_ld_clumping_population_parser(), "population",
    ) == argparse.SUPPRESS
    assert action_default(
        get_pipeline_genome_build_parser(), "genome_build",
    ) == argparse.SUPPRESS
    shared_build_factory = "postgwas.cli.common:get_pipeline_genome_build_parser"
    for module in ("annot_ldblock", "ld_clump", "qc_summary", "mixer"):
        assert shared_build_factory in REGISTRY.get(module).parser_factories

    imputation = configuration.modules.imputation
    imputation_parser = get_common_imputation_parser()
    pred_ld = imputation.engines.pred_ld
    assert action_default(imputation_parser, "imputation_engine") == imputation.engine
    assert action_default(imputation_parser, "imputation_r2_threshold") == (
        pred_ld.minimum_r2
    )
    assert action_default(imputation_parser, "imputation_minimum_maf") == (
        pred_ld.minimum_maf
    )
    assert action_default(imputation_parser, "ref") == pred_ld.mode
    assert action_default(imputation_parser, "corr_method") == (
        pred_ld.correlation_method
    )
    assert imputation_parser.parse_args(
        ["--correlation-method", "spearman"],
    ).corr_method == "spearman"
    assert action_default(
        get_imputation_population_parser(), "population",
    ) == imputation.population.value

    manhattan = configuration.modules.manhattan
    manhattan_parser = get_assoc_plot_parser()
    for destination, value in {
        "nauto": manhattan.autosome_count,
        "min_af": manhattan.minimum_af,
        "min_lp": manhattan.minimum_neglog10_p,
        "loglog_pval": manhattan.loglog_pvalue,
        "cyto_ratio": manhattan.cytoband_ratio,
        "max_height": manhattan.maximum_height,
        "spacing": manhattan.chromosome_spacing,
        "width": manhattan.width,
        "height": manhattan.height,
        "fontsize": manhattan.font_size,
        "allelic_shift": manhattan.allelic_shift,
        "csq": manhattan.flag_coding,
    }.items():
        assert action_default(manhattan_parser, destination) == value

    enrichment = configuration.modules.enrichment
    assert action_default(
        get_geneset_parser(add_help=True), "reference_set",
    ) == enrichment.reference_set


def test_pipeline_resolves_one_configured_genome_build_or_requires_override():
    configuration = load_configuration().model_copy(deep=True)
    args = argparse.Namespace()

    pipeline_cli._resolve_pipeline_genome_build(
        args, ("annot_ldblock", "ld_clump"), configuration,
    )
    assert args.genome_build == "GRCh37"

    configuration.modules.ld_clumping.genome_build = "GRCh38"
    with pytest.raises(
        ConfigurationError, match="incompatible genome builds",
    ) as captured:
        pipeline_cli._resolve_pipeline_genome_build(
            argparse.Namespace(),
            ("annot_ldblock", "ld_clump"),
            configuration,
        )
    assert "modules.ld_annotation.genome_build=GRCh37" in str(captured.value)
    assert "modules.ld_clumping.genome_build=GRCh38" in str(captured.value)

    explicit = argparse.Namespace(genome_build="GRCh38")
    pipeline_cli._resolve_pipeline_genome_build(
        explicit, ("annot_ldblock", "ld_clump"), configuration,
    )
    assert explicit.genome_build == "GRCh38"


def test_required_help_uses_one_red_prefix_at_the_beginning():
    markup = style_cli_requirement("Input file.", required=True)
    assert markup == "[bold red]Required:[/bold red] Input file."
    rendered = Text.from_markup(markup)
    assert rendered.plain == "Required: Input file."
    assert [str(span.style) for span in rendered.spans] == ["bold red"]
    assert style_cli_requirement(markup, required=True) == markup
    assert style_cli_requirement("REQUIRED. Input file.") == markup
    assert style_cli_requirement(
        "[bold bright_red]Required[/bold bright_red]: Input file."
    ) == markup


def test_required_help_marker_does_not_change_argparse_validation():
    parser = argparse.ArgumentParser(formatter_class=AlignedRichHelpFormatter)
    parser.add_argument("--from-config", help="Input that may be supplied by YAML.")
    mark_cli_required_help(parser, ("from_config",))

    assert parser.parse_args([]).from_config is None
    assert "Required: Input that may be supplied by YAML." in parser.format_help()


def test_conditional_requirement_help_does_not_change_parser_semantics():
    parser = argparse.ArgumentParser(formatter_class=AlignedRichHelpFormatter)
    parser.add_argument(
        "--conditional-input",
        default=argparse.SUPPRESS,
        help=help_with_conditional_requirement(
            "Mode-specific input.",
            "for method_x",
        ),
    )

    assert vars(parser.parse_args([])) == {}
    normalized = " ".join(parser.format_help().split())
    assert (
        "--conditional-input CONDITIONAL_INPUT Required: for method_x. "
        "Mode-specific input."
    ) in normalized


def test_public_help_excludes_applicability_annotations():
    invocations = list(PUBLIC_HELP_COMMANDS.values())
    invocations.extend(
        ("pipeline", "--modules", target, "--help")
        for target in PIPELINE_TARGETS
    )
    for arguments in invocations:
        normalized = " ".join(_help(arguments, 170).split()).lower()
        assert "applies to" not in normalized, arguments


def test_direct_help_keeps_conditional_requirement_annotations():
    expected = {
        "gcta_cojo": "--condition-snps PATH Required: for cond mode.",
        "flames": (
            "--cadd-file PATH Required: when --flames-cadd-mode local is selected."
        ),
    }
    for command, fragment in expected.items():
        normalized = " ".join(_help((command, "--help"), 170).split())
        assert fragment in normalized


def test_direct_help_excludes_pipeline_run_examples_and_generated_input_language():
    for command in ("kpops", "caldera"):
        help_text = _help((command, "--help"), 150)
        assert "postgwas pipeline" not in help_text

    for command in (
        "gcta_gene", "gcta_cojo", "magma", "magmacovar", "single_cell",
        "mixer",
    ):
        normalized = " ".join(_help((command, "--help"), 170).split()).lower()
        assert "pipeline mode supplies this automatically" not in normalized
        assert "pipeline supplies this automatically" not in normalized
        assert "pipeline run receives this file automatically" not in normalized


def test_pipeline_help_marks_registry_requirements_without_cli_duplication():
    normalized_help = " ".join(
        _help(("pipeline", "--modules", "heritability", "--help"), 160).split()
    )
    for option in (
        "--vcf PATH",
        "--dataset-id NAME",
        "--output-directory PATH",
        "--merge-alleles PATH",
        "--ref-ld-chr PATH",
        "--w-ld-chr PATH",
    ):
        assert "%s Required:" % option in normalized_help


def test_ldsc_prevalence_help_distinguishes_direct_and_pipeline_modes():
    direct_help = " ".join(_help(("heritability", "--help"), 160).split())
    pipeline_help = " ".join(
        _help(("pipeline", "--modules", "heritability", "--help"), 160).split()
    )

    assert "Direct mode requires it from --samp-prev or YAML" in direct_help
    assert "PostGWAS does not calculate or infer it" in direct_help
    assert "formatter-returned value" not in direct_help

    assert "Optional sample-prevalence override" in pipeline_help
    assert "PostGWAS uses the formatter-returned value" in pipeline_help
    assert "it does not recalculate it" in pipeline_help

    for help_text in (direct_help, pipeline_help):
        assert help_text.count("Default: unset") >= 2


def test_standalone_manual_requirement_is_visibly_marked():
    normalized_help = " ".join(_help(("imputation", "--help"), 160).split())
    assert "--pred-ld-input-directory PATH Required:" in normalized_help


def test_pipeline_error_renders_bracketed_external_tool_paths_literally(monkeypatch):
    stream = StringIO()
    monkeypatch.setattr(
        pipeline_cli,
        "console",
        Console(file=stream, force_terminal=False, color_system=None),
    )

    pipeline_cli._print_error(
        "Pipeline failed:",
        RuntimeError("GCTA failed for [/analysis/pathway.set]"),
    )

    assert stream.getvalue().strip() == (
        "Pipeline failed: GCTA failed for [/analysis/pathway.set]"
    )
