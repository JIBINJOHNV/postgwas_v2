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
    get_gcta_gene_pipeline_examples,
)
from postgwas.modules.harmonisation.cli import get_harmonisation_parser
from postgwas.modules.harmonisation.concordance.cli import get_validation_parser
from postgwas.modules.imputation.cli import build_parser as build_imputation_parser
from postgwas.modules.ld_annotation.cli import (
    build_parser as build_ld_annotation_parser,
)
from postgwas.modules.ld_clumping.cli import build_parser as build_ld_clumping_parser
from postgwas.modules.ldsc.cli import (
    build_parser as build_ldsc_parser,
    get_ldsc_pipeline_examples,
)
from postgwas.modules.magma.cli import build_parser as build_magma_parser
from postgwas.modules.magmacovar.cli import (
    build_parser as build_magmacovar_parser,
    get_magmacovar_pipeline_examples,
)
from postgwas.modules.manhattan.cli import build_parser as build_manhattan_parser
from postgwas.modules.mixer.cli import (
    build_parser as build_mixer_parser,
    get_mixer_parser,
)
from postgwas.modules.pops.cli import build_parser as build_pops_parser
from postgwas.modules.qc_summary.cli import build_parser as build_qc_parser
from postgwas.modules.single_cell.cli import build_parser as build_single_cell_parser
from postgwas.modules.fine_mapping.cli import (
    build_parser as build_finemap_parser,
    get_finemap_pipeline_examples,
)
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
    get_ld_clump_parser,
    get_ld_clumping_population_parser,
    get_pipeline_genome_build_parser,
    sumstat_summary_arg_parser,
)
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import (
    AlignedRichHelpFormatter,
    cli_option_value,
    cli_option_values,
    format_cli_choices,
    format_cli_default,
    format_cli_examples,
    format_cli_help_block,
    help_with_choices,
    help_with_conditional_requirement,
    help_with_default,
    mark_cli_required_help,
    style_cli_defaults,
    style_cli_requirement,
)
from postgwas.pipeline.registry import REGISTRY
from postgwas.pipeline import cli as pipeline_cli
from preflight_support import pipeline_input_vcf_evidence
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
        parser = builder()
        _assert_all_finite_choices_are_visible(parser, command)
        option_strings = {
            option
            for action in parser._actions
            for option in action.option_strings
        }
        assert "--magma" not in option_strings, command
        assert "--bcftools" not in option_strings, command


@pytest.mark.parametrize("width", (80, 120))
def test_filtering_boolean_help_uses_presence_flags_in_every_mode(width):
    help_pages = (
        _help(("sumstat_filter", "--help"), width),
        _help(("pipeline", "--modules", "sumstat_filter", "--help"), width),
    )
    for help_text in help_pages:
        normalized_help = " ".join(help_text.split())
        for option in (
            "include-indels",
            "remove-palindromic",
            "remove-mhc",
            "write-soft-filter-vcf",
        ):
            assert "--%s" % option in help_text
            assert "--%s BOOL" % option not in help_text
            assert "--no-%s" % option not in help_text
        assert "--mhc-chrom CHROM" in help_text
        assert "--mhc-start POSITION" in help_text
        assert "--mhc-end POSITION" in help_text
        assert "Default: GRCh37=6, GRCh38=6" in normalized_help
        assert (
            "Default: GRCh37=28477797, GRCh38=28510120"
            in normalized_help
        )
        assert (
            "Default: GRCh37=33448354, GRCh38=33480577"
            in normalized_help
        )
        assert "otherwise use mhc_regions (Default: unset)" not in normalized_help


def test_qc_and_filtering_share_public_policy_option_destinations():
    qc_actions = {
        option: action.dest
        for action in build_qc_parser()._actions
        for option in action.option_strings
    }
    filtering_actions = {
        option: action.dest
        for action in build_filtering_parser()._actions
        for option in action.option_strings
    }

    for option in (
        "--minimum-neglog10-p",
        "--minimum-maf",
        "--reference-af-column",
        "--maximum-af-difference",
        "--minimum-info",
        "--maximum-info",
        "--missing-pvalue-action",
        "--missing-af-action",
        "--missing-info-action",
        "--include-indels",
        "--remove-palindromic",
        "--palindromic-af-lower",
        "--palindromic-af-upper",
        "--remove-mhc",
    ):
        assert qc_actions[option] == filtering_actions[option]


@pytest.mark.parametrize("width", (80, 120))
def test_qc_policy_help_is_available_in_direct_and_pipeline_modes(width):
    help_pages = (
        _help(("qc", "--help"), width),
        _help(("pipeline", "--modules", "qc_summary", "--help"), width),
    )
    for help_text in help_pages:
        normalized_help = " ".join(help_text.split())
        for option in (
            "--minimum-neglog10-p",
            "--minimum-maf",
            "--reference-af-column",
            "--maximum-af-difference",
            "--minimum-info",
            "--maximum-info",
            "--missing-pvalue-action",
            "--missing-af-action",
            "--missing-info-action",
            "--palindromic-af-lower",
            "--palindromic-af-upper",
            "--sample-size-reference-quantile",
            "--sample-size-minimum-fraction",
        ):
            assert option in help_text
        for option in ("include-indels", "remove-palindromic", "remove-mhc"):
            assert "--%s" % option in help_text
            assert "--%s BOOL" % option not in help_text
            assert "--no-%s" % option not in help_text
        assert "Default: indels are not included" in normalized_help
        assert "Default: palindromic variants are not removed" in normalized_help
        assert "Default: MHC variants are not removed" in normalized_help
        assert "--minimum-info VALUE" in help_text
        assert "Default: 0.7" in normalized_help
        assert "--maximum-info VALUE" in help_text
        assert "Default: 1.05" in normalized_help
        assert "--genome-build" not in help_text
        assert "--bcftools" not in help_text


def test_combined_filtering_and_qc_pipeline_deduplicates_shared_policy_options():
    help_text = _help(
        ("pipeline", "--modules", "sumstat_filter", "qc_summary", "--help"),
        120,
    )

    for option in (
        "--minimum-neglog10-p",
        "--minimum-maf",
        "--reference-af-column",
        "--maximum-af-difference",
        "--minimum-info",
        "--maximum-info",
        "--missing-pvalue-action",
        "--missing-af-action",
        "--missing-info-action",
        "--include-indels",
        "--remove-palindromic",
        "--remove-mhc",
    ):
        assert option in help_text


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
            assert "--show-screen" not in help_text, command
            assert "--hide-screen" in help_text, command
            assert "--resume" in help_text, command
            assert "--no-resume" not in help_text, command
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
            "Required:",
            "Required:",
        ),
    ),
)
def test_cli_help_uses_shared_description_column_for_same_and_wrapped_rows(
    width, arguments, build_marker, population_marker,
):
    help_lines = _help(arguments, width).splitlines()
    build_line = next(
        line
        for line in help_lines
        if line.lstrip().startswith("--genome-build BUILD")
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
        if line.lstrip().startswith("--gcta-reference-population CODE")
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
    assert parser.parse_args(["--hide-screen"]).show_screen is False
    assert parser.parse_args(["--resume"]).resume is True
    assert parser.parse_args(["--overwrite"]).overwrite is True
    with pytest.raises(SystemExit):
        parser.parse_args(["--show-screen"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--resume", "--overwrite"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--no-resume"])

    actions = {
        action.dest: action
        for action in parser._actions
        if action.dest in {"show_screen", "resume", "overwrite"}
    }
    assert "[bold green]Default:[/bold green] [green]false[/green]" in (
        actions["show_screen"].help or ""
    )
    assert "[bold green]Default:[/bold green] [green]true[/green]" in (
        actions["resume"].help or ""
    )
    assert "[bold green]Default:[/bold green] [green]false[/green]" in (
        actions["overwrite"].help or ""
    )
    assert (
        "Reuse validated completed steps and continue from the first "
        "incomplete step"
    ) in (actions["resume"].help or "")
    assert (
        "Start the analysis from the first step and replace files created "
        "by the previous run"
    ) in (actions["overwrite"].help or "")


def test_every_pipeline_target_builds_context_specific_help():
    for module in PIPELINE_TARGETS:
        help_text = _help(("pipeline", "--modules", module, "--help"), 120)
        expected_usage = (
            "Usage: postgwas pipeline --modules finemap "
            "--clumping-methods METHOD [METHOD ...] "
            "--finemap-method METHOD [--help]"
            if module == "finemap"
            else "Usage: postgwas pipeline [options] --modules %s" % module
        )
        assert expected_usage in help_text
        assert "Workflow examples" in help_text
        assert "--show-screen" not in help_text
        assert "--hide-screen" in help_text
        assert "--resume" in help_text
        assert "--no-resume" not in help_text
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
        option_strings = {
            option
            for action in parser._actions
            for option in action.option_strings
        }
        assert "--magma" not in option_strings, module
        assert "--bcftools" not in option_strings, module
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

    assert "--resume" in help_text
    assert "--no-resume" not in help_text
    assert "--overwrite" in help_text
    assert "--resolve-variants-to-reference" in help_text
    assert "--magma-memory-per-worker-gb GB" in normalized_help
    assert "Default: 16 GB" in normalized_help
    assert "--duplicate-policy {err,lowest_p,remove}" in normalized_help
    assert "Default: lowest_p" in normalized_help
    assert "Available options: err, lowest_p, remove" in normalized_help
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
        lambda args, plan, configuration, **kwargs: executed.append(
            (args, plan, configuration, kwargs["preflight_evidence"])
        ),
    )
    validated_vcf = pipeline_input_vcf_evidence()["input_vcf"]
    monkeypatch.setattr(
        pipeline_cli,
        "_validate_pipeline_entry_vcf",
        lambda _args, _configuration: validated_vcf,
    )
    monkeypatch.setattr(
        pipeline_cli,
        "_run_pipeline_preflights",
        lambda _args, _modules, **_kwargs: {"input_vcf": validated_vcf},
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
    assert "--bcftools" not in help_text
    assert "--gcta PATH" not in help_text
    assert "External software:" not in help_text
    assert "GCTA binary:" not in help_text
    assert "Input files:" in help_text
    assert "Reference files:" in help_text
    assert "GCTA reference settings:" in help_text
    assert "GCTA analysis settings:" in help_text
    assert "GCTA variant matching:" in help_text
    assert "GCTA pathway conversion:" in help_text
    assert "GCTA result reporting:" in help_text
    section_offsets = [
        help_text.index(title)
        for title in (
            "Output:",
            "Input files:",
            "Reference files:",
            "GCTA reference settings:",
            "Variant identifiers:",
            "GCTA analysis settings:",
            "GCTA variant matching:",
            "GCTA pathway conversion:",
            "GCTA result reporting:",
        )
    ]
    assert section_offsets == sorted(section_offsets)
    assert "GCTA inputs:" not in help_text
    assert "--method METHOD" in help_text
    assert (
        "Available options: fastbat_gene, fastbat_segment, fastbat_set, mbat_combo"
        in normalized_help
    )
    assert "Default: unset" not in help_text
    for invocation in (
        "--gcta-reference-prefix PREFIX",
        "--genome-build BUILD",
        "--gcta-reference-population CODE",
    ):
        assert "%s Required:" % invocation in normalized_help
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("gcta_gene").required_options
    }
    assert required["gcta_reference_prefix"] == (
        "modules.gcta_gene.reference.prefix"
    )
    assert required["genome_build"] == "modules.gcta_gene.genome_build"
    assert required["gcta_reference_population"] == (
        "modules.gcta_gene.reference.population"
    )
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
    expected_examples = {
        "Run gene-based fastBAT:": "--method fastbat_gene",
        "Run fixed-segment fastBAT:": "--method fastbat_segment",
        "Run custom-set fastBAT from a prepared set list:": (
            "--fastbat-set-list pathways.fastbat.set"
        ),
        "Convert a GMT pathway file and run custom-set fastBAT:": (
            "--gmt pathways.gmt"
        ),
        "Run gene-based mBAT-combo:": "--method mbat_combo",
    }
    for label, argument in expected_examples.items():
        assert label in help_text
        assert argument in help_text


@pytest.mark.parametrize("width", (80, 120))
@pytest.mark.parametrize(
    ("method", "visible_options", "hidden_options", "example_labels"),
    (
        (
            "fastbat_gene",
            ("--gene-list", "--gene-window-kb"),
            (
                "--fastbat-set-list", "--gmt", "--fastbat-segment-size-kb",
                "--mbat-svd-gamma",
                "--print-component-p-values",
            ),
            ("Run gene-based fastBAT:",),
        ),
        (
            "fastbat_segment",
            ("--fastbat-segment-size-kb",),
            (
                "--gene-list", "--gene-window-kb", "--fastbat-set-list",
                "--gmt", "--mbat-svd-gamma",
                "--print-component-p-values",
            ),
            ("Run fixed-segment fastBAT:",),
        ),
        (
            "fastbat_set",
            (
                "--gene-list", "--gene-window-kb", "--fastbat-set-list",
                "--gmt", "--gmt-empty-pathway-policy",
                "--fastbat-oversized-set-policy",
            ),
            (
                "--fastbat-segment-size-kb", "--mbat-svd-gamma",
                "--print-component-p-values",
            ),
            (
                "Run custom-set fastBAT from a prepared set list:",
                "Convert a GMT pathway file and run custom-set fastBAT:",
            ),
        ),
        (
            "mbat_combo",
            (
                "--gene-list", "--gene-window-kb", "--mbat-svd-gamma",
                "--frequency-difference-max", "--print-component-p-values",
            ),
            ("--fastbat-set-list", "--gmt", "--fastbat-segment-size-kb"),
            ("Run gene-based mBAT-combo:",),
        ),
    ),
)
def test_pipeline_gcta_help_is_filtered_to_the_selected_method(
    width, method, visible_options, hidden_options, example_labels,
):
    help_text = _help(
        ("pipeline", "--modules", "gcta_gene", "--method", method, "--help"),
        width,
    )
    normalized_help = " ".join(help_text.split())

    assert (
        "This help page is filtered for the selected procedure:"
        in normalized_help
    )
    assert "%s: " % method in help_text
    assert "--method METHOD" in help_text
    assert (
        "Available options: fastbat_gene, fastbat_segment, fastbat_set, mbat_combo"
        in normalized_help
    )
    for option in (
        "--vcf", "--gcta-reference-prefix", "--genome-build",
        "--gcta-reference-population", "--gcta-reference-maf-min",
        "--fastbat-ld-cutoff", "--write-snpset", "--gcta-top-results",
        "--frequency-difference-max",
    ):
        assert option in help_text
    for option in visible_options:
        assert option in help_text
    for option in hidden_options:
        assert option not in help_text

    all_example_labels = {
        "Run gene-based fastBAT:",
        "Run fixed-segment fastBAT:",
        "Run custom-set fastBAT from a prepared set list:",
        "Convert a GMT pathway file and run custom-set fastBAT:",
        "Run gene-based mBAT-combo:",
    }
    for label in example_labels:
        assert label in help_text
    for label in all_example_labels - set(example_labels):
        assert label not in help_text


@pytest.mark.parametrize(
    ("source_arguments", "visible_options", "hidden_options", "example_label"),
    (
        (
            ("--fastbat-set-list", "pathways.fastbat.set"),
            (
                "--fastbat-set-list", "--gmt-empty-pathway-policy",
                "--fastbat-oversized-set-policy",
            ),
            (
                "--gmt", "--gene-list", "--gene-window-kb",
                "--gmt-chromosome-label-policy", "--minimum-gene-id-overlap",
            ),
            "Run custom-set fastBAT from a prepared set list:",
        ),
        (
            ("--gmt", "pathways.gmt"),
            (
                "--gmt", "--gene-list", "--gene-window-kb",
                "--gmt-chromosome-label-policy", "--minimum-gene-id-overlap",
                "--gmt-empty-pathway-policy", "--fastbat-oversized-set-policy",
            ),
            ("--fastbat-set-list",),
            "Convert a GMT pathway file and run custom-set fastBAT:",
        ),
    ),
)
def test_pipeline_fastbat_set_help_is_filtered_to_the_selected_input_source(
    source_arguments, visible_options, hidden_options, example_label,
):
    help_text = _help(
        (
            "pipeline", "--modules", "gcta_gene", "--method", "fastbat_set",
            *source_arguments, "--help",
        ),
        120,
    )
    normalized_help = " ".join(help_text.split())
    displayed_options = {
        match.group(1)
        for line in help_text.splitlines()
        for match in [re.match(r"^\s+(--[a-z0-9-]+)", line)]
        if match is not None
    }

    for option in visible_options:
        assert option in displayed_options
    for option in hidden_options:
        assert option not in displayed_options
    selected_option = source_arguments[0]
    selected_metavar = "PATH"
    assert "%s %s Required:" % (selected_option, selected_metavar) in normalized_help
    if selected_option == "--gmt":
        assert "--gene-list PATH Required:" in normalized_help
    assert help_text.count("Run custom-set fastBAT from a prepared set list:") == (
        1 if example_label.startswith("Run custom-set") else 0
    )
    assert help_text.count(
        "Convert a GMT pathway file and run custom-set fastBAT:"
    ) == (1 if example_label.startswith("Convert") else 0)


def test_pipeline_gcta_missing_external_inputs_stop_before_execution(tmp_path):
    input_vcf = Path("tests/fixtures/formatting/harmonised.vcf").resolve()
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "postgwas",
            "pipeline",
            "--modules",
            "gcta_gene",
            "--vcf",
            str(input_vcf),
            "--dataset-id",
            "STUDY",
            "--output-directory",
            str(tmp_path / "results"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )

    assert completed.returncode == 2
    terminal = completed.stdout + completed.stderr
    expectations = {
        "--gcta-reference-prefix": "modules.gcta_gene.reference.prefix",
        "--genome-build": "modules.gcta_gene.genome_build",
        "--gcta-reference-population": (
            "modules.gcta_gene.reference.population"
        ),
    }
    terminal = " ".join(terminal.split())
    screen_log = " ".join((
        tmp_path / "results" / "run_metadata" / "screen.log"
    ).read_text(encoding="utf-8").split())
    for option, configuration_path in expectations.items():
        message = (
            "Required argument not provided: %s. Provide %s VALUE or set %s "
            "in the run configuration."
            % (option, option, configuration_path)
        )
        assert message in terminal
        assert message in screen_log
    generated_files = {
        path.relative_to(tmp_path / "results").as_posix()
        for path in (tmp_path / "results").rglob("*")
        if path.is_file()
    }
    assert generated_files == {"run_metadata/screen.log"}


def test_gcta_pipeline_examples_use_public_options_and_mode_specific_inputs():
    option_strings = {
        option
        for action in build_gcta_gene_parser()._actions
        for option in action.option_strings
    }
    pipeline_only_options = {"--modules", "--vcf"}
    examples = get_gcta_gene_pipeline_examples()

    assert len(examples) == 5
    for _label, command, arguments in examples:
        assert command == "postgwas pipeline"
        example_options = {argument.split()[0] for argument in arguments}
        assert example_options <= option_strings | pipeline_only_options
        assert {
            "--modules",
            "--vcf",
            "--gcta-reference-prefix",
            "--genome-build",
            "--gcta-reference-population",
            "--dataset-id",
            "--output-directory",
        } <= example_options

        method = next(
            argument.split()[1]
            for argument in arguments
            if argument.startswith("--method ")
        )
        if method in {"fastbat_gene", "mbat_combo"}:
            assert "--gene-list" in example_options
        elif method == "fastbat_segment":
            assert "--gene-list" not in example_options
        elif method == "fastbat_set":
            set_sources = example_options & {"--fastbat-set-list", "--gmt"}
            assert len(set_sources) == 1
            assert ("--gene-list" in example_options) == ("--gmt" in set_sources)


def test_pipeline_gcta_cojo_help_marks_external_inputs_and_all_four_modes():
    help_text = _help(("pipeline", "--modules", "gcta_cojo", "--help"), 120)
    normalized = " ".join(help_text.split())

    assert "--gcta-cojo-input-file" not in help_text
    for invocation in (
        "--vcf PATH",
        "--cojo-reference-prefix PREFIX",
        "--genome-build BUILD",
        "--cojo-reference-population CODE",
    ):
        assert "%s Required:" % invocation in normalized
    for label in (
        "Select independent signals when no conditioning SNPs are supplied:",
        "Select up to a requested number of independent signals:",
        "Estimate a specified SNP set jointly:",
        "Run conditional association on an explicit SNP list:",
    ):
        assert label in help_text
    assert "adjusted P-value" not in help_text
    assert "adjusted-p" not in help_text


def test_pipeline_pops_help_shows_required_external_resources():
    direct_help = build_pops_parser().format_help()
    help_text = _help(("pipeline", "--modules", "pops", "--help"), 120)
    normalized_help = " ".join(help_text.split())

    assert "PoPS does not use MAGMACOVAR .gsa.out results" in direct_help
    assert "technical covariates read from MAGMA .genes.raw" in direct_help
    assert "--gene-universe-policy POLICY" in direct_help
    assert "--gene-universe-policy intersect" in direct_help
    assert "Run MAGMA gene analysis followed by PoPS prioritisation:" in help_text
    assert (
        "Run MAGMA and explicitly intersect incompatible gene-score genes:"
        in help_text
    )
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
    for invocation in (
        "--feature-matrix-prefix PREFIX",
        "--pops-gene-location-file PATH",
        "--pops-genome-build BUILD",
    ):
        assert "%s Required:" % invocation in normalized_help
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("pops").required_options
    }
    assert required["feature_matrix_prefix"] == (
        "modules.pops.feature_matrix_prefix"
    )
    assert required["pops_gene_location_file"] == (
        "modules.pops.gene_location_file"
    )
    assert required["pops_genome_build"] == "modules.pops.genome_build"
    assert "--run-config analysis_pipeline.yaml" not in help_text


def test_pops_help_groups_outputs_inputs_and_references_before_settings(
    monkeypatch,
):
    direct_parser = build_pops_parser()
    captured = {}

    def capture_help(modules, parser, *, error=False):
        captured["parser"] = parser

    monkeypatch.setattr(pipeline_cli, "print_full_pipeline_help", capture_help)
    monkeypatch.setattr(
        sys,
        "argv",
        ["postgwas pipeline", "--modules", "pops", "--help"],
    )
    with pytest.raises(SystemExit) as stopped:
        pipeline_cli.main()
    assert stopped.value.code == 0

    for mode, parser in (
        ("direct", direct_parser),
        ("pipeline", captured["parser"]),
    ):
        visible_groups = [
            group
            for group in parser._action_groups
            if group.description or any(
                action.help is not argparse.SUPPRESS
                for action in group._group_actions
            )
        ]
        assert [group.title for group in visible_groups[:4]] == [
            "options",
            "Output",
            "Input files",
            "Reference files",
        ], mode
        destinations = {
            group.title: [
                action.dest
                for action in group._group_actions
                if action.help is not argparse.SUPPRESS
            ]
            for group in visible_groups
        }
        assert destinations["Output"] == [
            "dataset_id",
            "output_directory",
            "save_matrix_files",
            "save_matrix_files",
        ], mode
        expected_inputs = [
            "target_score_file",
            "target_covariates_file",
            "target_error_covariance_file",
        ]
        if mode == "direct":
            expected_inputs[0:0] = [
                "magma_association_prefix",
                "magma_annotated_results_file",
            ]
        else:
            expected_inputs = ["vcf"]
            pipeline_help = parser.format_help()
            for option in (
                "--target-score-file",
                "--target-covariates-file",
                "--target-error-covariance-file",
            ):
                assert option not in pipeline_help
        assert destinations["Input files"] == expected_inputs, mode
        assert "PoPS input" not in destinations, mode
        assert "PoPS input compatibility" not in destinations, mode
        if mode == "direct":
            input_group = next(
                group
                for group in visible_groups
                if group.title == "Input files"
            )
            assert "requires exactly one gene-score source" in (
                input_group.description or ""
            )
            direct_help = parser.format_help()
            compact_direct_help = " ".join(direct_help.split())
            assert "--target-score-file" in direct_help
            assert "Custom gene-score table" in compact_direct_help
            assert (
                "The option name is retained for compatibility"
                in compact_direct_help
            )

        expected_references = [
            "feature_matrix_prefix",
            "pops_gene_location_file",
            "feature_subset_file",
            "control_features_file",
        ]
        if mode == "pipeline":
            expected_references[:0] = [
                "magma_ld_reference",
                "gene_location_file",
                "gene_set_file",
            ]
        assert destinations["Reference files"] == expected_references, mode
        for title in (
            "PoPS covariate settings",
            "PoPS feature-selection settings",
            "PoPS training settings",
            "PoPS model settings",
        ):
            assert visible_groups.index(
                next(group for group in visible_groups if group.title == title)
            ) > 3, (mode, title)
        ordered_setting_titles = [
            title
            for title in (
                "PoPS resource settings",
                "Variant identifiers",
                "MAGMA analysis settings",
                "MAGMA analysis scope",
                "PoPS covariate settings",
                "PoPS feature-selection settings",
                "PoPS training settings",
                "PoPS model settings",
                "Performance",
                "Screen output",
                "Run continuation",
                "Configuration",
                "Choose analyses",
            )
            if title in destinations
        ]
        assert [
            group.title
            for group in visible_groups[4:]
            if group.title in ordered_setting_titles
        ] == ordered_setting_titles, mode


def test_finemap_pipeline_help_routes_through_method_selection():
    help_text = _help(("pipeline", "--modules", "finemap", "--help"), 120)
    normalized = " ".join(help_text.split())

    assert "Fine-mapping supports multiple workflow combinations." in help_text
    assert "--clumping-methods METHOD [METHOD ...]" in help_text
    assert "--finemap-method METHOD" in help_text
    assert "Available options: region, standard, cojo-slct" in normalized
    assert "Available options: susie, finemap" in normalized
    assert "Default: region standard" in normalized
    assert "Default: susie" in normalized
    assert "--run-config PATH" in help_text
    assert "--resume" in help_text
    assert "--overwrite" in help_text
    assert "--hide-screen" in help_text
    assert "Choose the fine-mapping workflow" in help_text
    assert "Steps that will run" not in help_text
    for hidden in (
        "--vcf PATH",
        "--ld-folder PATH",
        "--cojo-reference-prefix PREFIX",
        "--L N",
        "--n-causal-snps N",
    ):
        assert hidden not in help_text
    assert (
        "Show every option for standard clumping with SuSiE-RSS:"
        in help_text
    )
    assert (
        "Show every option for standard clumping with FINEMAP:"
        in help_text
    )


@pytest.mark.parametrize("width", (80, 120))
@pytest.mark.parametrize(
    (
        "engine", "visible_options", "hidden_options", "plink_default",
    ),
    (
        (
            "susie",
            ("--L N",),
            (
                "--n-causal-snps N",
                "--ldstore-timeout-seconds SECONDS",
                "--sss",
                "--cond",
            ),
            "plink",
        ),
        (
            "finemap",
            ("--n-causal-snps N", "--sss", "--cond"),
            ("--L N", "--susie-timeout-seconds SECONDS"),
            "plink2",
        ),
    ),
)
def test_finemap_pipeline_help_is_complete_for_selected_engine(
    width, engine, visible_options, hidden_options, plink_default,
):
    help_text = _help(
        (
            "pipeline", "--modules", "finemap",
            "--clumping-methods", "standard",
            "--finemap-method", engine,
            "--help",
        ),
        width,
    )
    normalized = " ".join(help_text.split())

    assert "Selected fine-mapping workflow" in help_text
    assert "LD-clumping methods : standard r²" in help_text
    assert "Fine-mapping method : %s" % (
        "SuSiE-RSS" if engine == "susie" else "FINEMAP"
    ) in help_text
    assert "1) Standard LD clumping (r²)" in help_text
    assert "2) %s input preparation" % (
        "SuSiE-RSS" if engine == "susie" else "FINEMAP"
    ) in help_text
    assert "3) %s fine-mapping" % (
        "SuSiE-RSS" if engine == "susie" else "FINEMAP"
    ) in help_text
    for option in (
        "--vcf PATH",
        "--ld-folder PATH",
        "--ld-window-kb KB",
        "--window-kb INT",
        "--finemap-ld-reference PREFIX",
        "--plink PATH",
        *visible_options,
    ):
        assert option in help_text
    for option in (
        *hidden_options,
        "--clumping-methods METHOD",
        "--finemap-method METHOD",
        "--ld-region-dir PATH",
        "--ld-block-populations POPULATION",
    ):
        assert option not in help_text
    for cojo_only in (
        "--cojo-reference-prefix PREFIX",
        "--cojo-p P",
        "--cojo-merge-dist BP",
    ):
        assert cojo_only not in help_text
    assert "Default: %s" % plink_default in normalized
    # The resolved values remain in usage/examples, while their selector rows
    # are removed from the already selected option list.
    assert "--clumping-methods standard" in help_text
    assert "--finemap-method %s" % engine in help_text


def test_finemap_pipeline_help_includes_cojo_only_when_selected():
    help_text = _help(
        (
            "pipeline", "--modules", "finemap",
            "--clumping-methods", "standard", "cojo-slct",
            "--finemap-method", "susie", "--help",
        ),
        120,
    )

    for option in (
        "--cojo-reference-prefix PREFIX",
        "--cojo-p P",
        "--cojo-merge-dist BP",
    ):
        assert option in help_text
    assert "1) SuSiE-RSS input preparation" in help_text
    assert "2) LD clumping (standard r² + GCTA-COJO)" in help_text
    assert "3) SuSiE-RSS fine-mapping" in help_text
    assert "1) annot_ldblock" not in help_text


def test_finemap_pipeline_help_adds_ld_annotation_only_for_region_clumping():
    help_text = _help(
        (
            "pipeline", "--modules", "finemap",
            "--clumping-methods", "region", "standard",
            "--finemap-method", "susie", "--help",
        ),
        120,
    )

    assert "1) annot_ldblock" in help_text
    assert "2) LD clumping (annotated-region + standard r²)" in help_text
    assert "3) SuSiE-RSS input preparation" in help_text
    assert "4) SuSiE-RSS fine-mapping" in help_text
    assert "--ld-region-dir PATH" in help_text
    assert "--ld-block-populations POPULATION" in help_text


def test_finemap_pipeline_bare_selection_stops_before_analysis():
    completed = subprocess.run(
        [sys.executable, "-m", "postgwas", "pipeline", "--modules", "finemap"],
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )

    assert completed.returncode == 2
    terminal = completed.stdout + completed.stderr
    assert "Fine-mapping workflow selection is required." in terminal
    assert "--clumping-methods METHOD [METHOD ...]" in terminal
    assert "--finemap-method METHOD" in terminal
    assert "No analysis was started." in terminal
    assert "Starting Execution Chain" not in terminal


def test_selected_finemap_pipeline_missing_inputs_gets_filtered_help():
    completed = subprocess.run(
        [
            sys.executable, "-m", "postgwas", "pipeline",
            "--modules", "finemap",
            "--clumping-methods", "standard",
            "--finemap-method", "susie",
        ],
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )

    assert completed.returncode == 2
    terminal = completed.stdout + completed.stderr
    assert (
        "Required argument not provided: --vcf. Provide --vcf VALUE."
        in " ".join(terminal.split())
    )
    for option, configuration_path in (
        ("--ld-folder", "modules.ld_clumping.reference.directory"),
        (
            "--finemap-ld-reference",
            "modules.fine_mapping.input.ld_reference_prefix",
        ),
    ):
        assert (
            "Required argument not provided: %s. Provide %s VALUE or set %s "
            "in the run configuration."
            % (option, option, configuration_path)
        ) in " ".join(terminal.split())
    assert "--ld-region-dir" not in terminal
    assert "--L N" in terminal
    assert "--n-causal-snps N" not in terminal
    assert "Starting Execution Chain" not in terminal


def test_finemap_pipeline_requires_standard_as_the_locus_source():
    completed = subprocess.run(
        [
            sys.executable, "-m", "postgwas", "pipeline",
            "--modules", "finemap",
            "--clumping-methods", "region",
            "--finemap-method", "susie",
            "--help",
        ],
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )

    assert completed.returncode == 2
    terminal = completed.stdout + completed.stderr
    assert "Fine-mapping requires standard LD clumping" in " ".join(
        terminal.split()
    )
    assert "No analysis was started." not in terminal


def test_finemap_pipeline_help_resolves_methods_from_run_configuration(tmp_path):
    run_config = tmp_path / "finemap_help.yaml"
    run_config.write_text(
        "modules:\n"
        "  ld_clumping:\n"
        "    methods: [standard]\n"
        "  fine_mapping:\n"
        "    engine: finemap\n",
        encoding="utf-8",
    )

    help_text = _help(
        (
            "pipeline", "--modules", "finemap",
            "--run-config", str(run_config), "--help",
        ),
        120,
    )

    assert "LD-clumping methods : standard r²" in help_text
    assert "Fine-mapping method : FINEMAP" in help_text
    assert "--n-causal-snps N" in help_text
    assert "--L N" not in help_text


def test_context_help_option_reader_supports_scalar_and_multivalue_forms():
    arguments = (
        "--clumping-methods=standard", "cojo-slct",
        "--finemap-method", "susie", "--help",
    )

    assert cli_option_values(arguments, "--clumping-methods") == (
        "standard", "cojo-slct",
    )
    assert cli_option_value(arguments, "--finemap-method") == "susie"


def test_finemap_help_has_runnable_direct_and_selected_pipeline_examples():
    direct_help = _help(("finemap", "--help"), 120)
    pipeline_help = _help(
        (
            "pipeline", "--modules", "finemap",
            "--clumping-methods", "standard",
            "--finemap-method", "susie", "--help",
        ),
        120,
    )

    assert "--window-kb INT" in direct_help
    assert "--window-kb INT" in pipeline_help
    assert "--ld-window-kb KB" in pipeline_help
    assert "Default: 500" in direct_help
    assert "Default: 500" in pipeline_help
    assert "--locus-file" not in pipeline_help

    for option in (
        "--susie-input-file formatted/STUDY_susie.tsv.gz",
        "--locus-file loci.tsv",
        "--locus-type range",
        "--window-kb 0",
        "--finemap-ld-reference reference/1000G_EUR",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in direct_help
    for option in (
        "--vcf study_GRCh37.vcf.gz",
        "--clumping-methods standard",
        "--genome-build GRCh37",
        "--ld-folder reference/pairwise_ld",
        "--population EUR",
        "--variant-id-type unique",
        "--finemap-method susie",
        "--finemap-ld-reference reference/1000G_EUR",
        "--plink plink",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--ld-region-dir reference/ld_blocks" not in pipeline_help
    assert "--ld-block-populations EUR" not in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("finemap").required_options
    }
    assert required["finemap_ld_reference"] == (
        "modules.fine_mapping.input.ld_reference_prefix"
    )
    assert required["ld_folder"] == (
        "modules.ld_clumping.reference.directory"
    )

    for engine in ("susie", "finemap"):
        examples = get_finemap_pipeline_examples(
            clumping_methods=("standard",), engine=engine,
        )
        assert len(examples) == 1
        _label, command, arguments = examples[0]
        assert command == "postgwas pipeline"
        assert "--clumping-methods standard" in arguments
        assert "--finemap-method %s" % engine in arguments


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
            "--gene-location-file reference/PoPS_GRCh37_strand_aware.loc"
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
    assert "--kpops-gene-universe-policy POLICY" in pipeline_help
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


def test_kpops_and_caldera_help_group_outputs_inputs_and_references_before_settings(
    monkeypatch,
):
    def pipeline_parser(target):
        captured = {}

        def capture_help(modules, parser, *, error=False):
            captured["parser"] = parser

        monkeypatch.setattr(
            pipeline_cli,
            "print_full_pipeline_help",
            capture_help,
        )
        monkeypatch.setattr(
            sys,
            "argv",
            ["postgwas pipeline", "--modules", target, "--help"],
        )
        with pytest.raises(SystemExit) as stopped:
            pipeline_cli.main()
        assert stopped.value.code == 0
        return captured["parser"]

    cases = (
        (
            "kpops-direct",
            build_kpops_parser(),
            ["magma_association_prefix"],
            ["kpops_gene_annotation_file", "kernel_matrix_prefix"],
        ),
        (
            "kpops-pipeline",
            pipeline_parser("kpops"),
            ["vcf"],
            [
                "magma_ld_reference",
                "gene_location_file",
                "gene_set_file",
                "kpops_gene_annotation_file",
                "kernel_matrix_prefix",
            ],
        ),
        (
            "caldera-direct",
            build_caldera_parser(),
            ["pops_file", "credible_set_file"],
            ["caldera_repository"],
        ),
        (
            "caldera-pipeline",
            pipeline_parser("caldera"),
            ["vcf"],
            [
                "magma_ld_reference",
                "gene_location_file",
                "gene_set_file",
                "feature_matrix_prefix",
                "pops_gene_location_file",
                "feature_subset_file",
                "control_features_file",
                "ld_region_dir",
                "cojo_reference_prefix",
                "ld_folder",
                "finemap_ld_reference",
                "caldera_repository",
            ],
        ),
    )

    for mode, parser, expected_inputs, expected_references in cases:
        visible_groups = [
            group
            for group in parser._action_groups
            if group.description or any(
                action.help is not argparse.SUPPRESS
                for action in group._group_actions
            )
        ]
        assert [group.title for group in visible_groups[:4]] == [
            "options",
            "Output",
            "Input files",
            "Reference files",
        ], mode
        destinations = {
            group.title: [
                action.dest
                for action in group._group_actions
                if action.help is not argparse.SUPPRESS
            ]
            for group in visible_groups
        }
        assert destinations["Input files"] == expected_inputs, mode
        assert destinations["Reference files"] == expected_references, mode
        for legacy_title in (
            "Input file",
            "MAGMA inputs",
            "PoPS input",
            "PoPS input compatibility",
            "K-POPS inputs",
            "CALDERA inputs",
        ):
            assert legacy_title not in destinations, (mode, legacy_title)

        reference_index = next(
            index
            for index, group in enumerate(visible_groups)
            if group.title == "Reference files"
        )
        for group in visible_groups[reference_index + 1:]:
            assert group.title not in {"Output", "Input files", "Reference files"}

    caldera_pipeline_help = cases[-1][1].format_help()
    for generated_input in (
        "--pops-file",
        "--credible-set-file",
        "--target-score-file",
        "--target-covariates-file",
        "--target-error-covariance-file",
    ):
        assert generated_input not in caldera_pipeline_help


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
    missing_values_action = next(
        action
        for action in parser._actions
        if action.dest == "covariate_missing_values"
    )
    missing_genes_action = next(
        action
        for action in parser._actions
        if action.dest == "covariate_missing_genes"
    )
    maximum_missing_action = next(
        action
        for action in parser._actions
        if action.dest == "covariate_max_miss"
    )
    assert "Accepted modifiers: %s" % ", ".join(configured.model_modifiers) in (
        model_action.help or ""
    )
    assert tuple(direction_action.choices or ()) == (
        "two-sided",
        "greater",
        "smaller",
    )
    assert tuple(missing_values_action.choices or ()) == (
        "drop",
        "median",
        "mean",
    )
    assert tuple(missing_genes_action.choices or ()) == ("drop", "fill")
    for action in (
        missing_values_action,
        maximum_missing_action,
        missing_genes_action,
    ):
        assert action.default == argparse.SUPPRESS
    assert "--covariate-direction DIRECTION" in help_text
    for option in (
        "--covariate-missing-values ACTION",
        "--covariate-max-miss FRACTION",
        "--covariate-missing-genes ACTION",
    ):
        assert option in help_text
    assert "Available options: drop, median, mean" in normalized_help
    assert "Available options: drop, fill" in normalized_help
    assert "Values must be between 0 and MAGMA's maximum 0.2" in normalized_help
    assert "Default: median" in normalized_help
    assert "Default: 0.05" in normalized_help
    assert "Default: drop" in normalized_help
    for required_option in (
        "--magma-gene-results-file PATH",
        "--covariates PATH",
        "--dataset-id NAME",
        "--output-directory PATH",
    ):
        assert "%s Required:" % required_option in normalized_help

    for placeholder, description in configured.model_placeholders.items():
        assert "<%s> %s" % (placeholder, description) in normalized_help
    for specification in configured.model_modifiers.values():
        for syntax in specification.syntax:
            assert syntax in normalized_help
        assert specification.description in normalized_help
        assert specification.published_example in normalized_help
        assert specification.example in normalized_help
    assert "How to choose a MAGMA model" in help_text
    assert "Pre-MAGMA safety checks" in help_text
    assert "Property missingness guard" in normalized_help
    assert "stops with the property name" in normalized_help
    assert "missing fraction > max-miss" in normalized_help
    assert "missing-genes=fill" in normalized_help
    assert "All eligible genes from the .genes.raw file" in normalized_help
    assert "missing-genes=drop" in normalized_help
    assert "Only gene IDs overlapping" in normalized_help
    assert "direction-covar=<DIRECTION>" in normalized_help
    assert "not the blanket direction=" in normalized_help
    assert "MAGMA uses smaller, not less" in normalized_help
    assert (
        "Default: max-miss=0.05, missing-genes=drop, missing-values=median"
        in normalized_help
    )
    assert "Question: Does each property associate" in normalized_help
    assert "Model: omit --covariate-model" in normalized_help
    assert "Direction: --covariate-direction two-sided" in normalized_help
    assert "Question: Is a tissue positively associated" in normalized_help
    assert "Model: --covariate-model condition-hide=Average" in normalized_help
    assert "Direction: --covariate-direction greater" in normalized_help
    assert "Advanced MAGMA model modifiers" in help_text
    assert "Example property names are illustrative" in normalized_help
    assert normalized_help.count("Published GWAS example") == len(
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


@pytest.mark.parametrize("width", (80, 120))
def test_magmacovar_pipeline_help_is_concise_gene_only_and_aligned(width):
    help_text = _help(("pipeline", "--modules", "magmacovar", "--help"), width)
    normalized_help = " ".join(help_text.split())

    assert "How to choose a MAGMA model" not in help_text
    assert "Pre-MAGMA safety checks" not in help_text
    assert "Property missingness guard" not in normalized_help
    assert "Advanced MAGMA model modifiers" not in help_text
    assert "Published GWAS example" not in help_text
    for modifier in load_configuration().modules.magmacovar.model_modifiers:
        assert "MAGMA model modifier: %s" % modifier not in help_text
    assert "--covariates PATH Required:" in normalized_help
    assert "--covariate-direction DIRECTION" in help_text
    assert "--covariate-missing-values ACTION" in help_text
    assert "--covariate-max-miss FRACTION" in help_text
    assert "--covariate-missing-genes ACTION" in help_text
    assert "Available options: two-sided, greater, smaller" in normalized_help
    assert "Run marginal gene-property tests, as in the original FLAMES README:" in (
        help_text
    )
    assert "Run FUMA-style tissue-specific gene-property tests:" in help_text
    assert "--magma-ld-reference reference/1000G_EUR" in help_text
    assert "--gene-location-file reference/NCBI37.3.gene.loc" in help_text
    for pathway_option in (
        "--gene-set-file",
        "--minimum-gene-id-overlap",
        "--gene-set-identifier-mismatch",
        "--alternate-gene-id-duplicate-policy",
        "--gene-location-alternate-id-type",
    ):
        assert pathway_option not in help_text
    assert "pathway" not in help_text.lower()
    assert "gene-set" not in help_text.lower()
    assert (
        "Run the MAGMA gene-association analysis required by MAGMAcovar."
        in help_text
    )
    assert "--covariates reference/GTEx/gtex_v8_ts_avg_log2TPM.txt" in help_text
    assert (
        "--covariates reference/GTEx/gtex_v8_ts_general_avg_log2TPM.txt"
        in help_text
    )
    assert "--covariate-model condition-hide=Average" in help_text
    assert "--covariate-direction greater" in help_text
    assert "--run-config analysis_pipeline.yaml" not in help_text

    lines = help_text.splitlines()
    covariates_line = next(line for line in lines if "--covariates PATH" in line)
    description_column = next(
        line.index("Required:")
        for line in lines[lines.index(covariates_line):lines.index(covariates_line) + 2]
        if "Required:" in line
    )
    model_index = next(
        index for index, line in enumerate(lines) if "--covariate-model MODEL" in line
    )
    model_description = next(
        line
        for line in lines[model_index:model_index + 2]
        if "Choose how MAGMA" in line
    )
    assert model_description.index("Choose how MAGMA") == description_column


def test_magmacovar_model_reference_uses_shared_help_colours():
    parser = build_magmacovar_parser()
    descriptions = "\n".join(
        group.description or ""
        for group in parser._action_groups
        if group.title in {
            "Pre-MAGMA safety checks",
            "How to choose a MAGMA model",
            "Advanced MAGMA model modifiers",
        } or group.title.startswith("MAGMA model modifier:")
    )
    rendered = Text.from_markup(descriptions)
    styles = {str(span.style) for span in rendered.spans}

    assert "bold cyan" in styles
    assert "cyan" in styles
    assert "bold bright_yellow" in styles
    assert "bold magenta" in styles
    assert "[bold" not in _help(("magmacovar", "--help"), 120)


def test_cli_examples_use_the_shared_high_contrast_palette():
    rendered = Text.from_markup(format_cli_examples(
        (
            "Run an analysis:",
            "postgwas pipeline",
            ("--modules magmacovar", "--vcf study.vcf.gz"),
        ),
        notes=("Replace the example paths.",),
        title="Workflow examples",
    ))
    styles = {str(span.style) for span in rendered.spans}

    assert {
        "bold cyan",
        "bold bright_yellow",
        "bold magenta",
        "cyan",
        "dim",
    } <= styles
    assert "[bold" not in rendered.plain


def test_magmacovar_pipeline_examples_use_the_public_parser_options():
    option_strings = {
        option
        for action in build_magmacovar_parser()._actions
        for option in action.option_strings
    }
    pipeline_only_options = {
        "--modules",
        "--vcf",
        "--magma-ld-reference",
        "--gene-location-file",
    }

    for _label, _command, arguments in get_magmacovar_pipeline_examples():
        example_options = {argument.split()[0] for argument in arguments}
        assert example_options <= option_strings | pipeline_only_options
        for required in (
            "--modules",
            "--vcf",
            "--magma-ld-reference",
            "--gene-location-file",
            "--covariates",
            "--dataset-id",
            "--output-directory",
        ):
            assert required in example_options


def test_explicit_magma_and_magmacovar_targets_retain_magma_pathway_help():
    help_text = _help(
        ("pipeline", "--modules", "magma", "magmacovar", "--help"),
        120,
    )

    assert "--gene-set-file PATH" in help_text
    assert "--minimum-gene-id-overlap FRACTION" in help_text
    assert "Run MAGMA gene and gene-set association analysis." in help_text
    assert (
        "Run the MAGMA gene-association analysis required by MAGMAcovar."
        not in help_text
    )


def test_magmacovar_pipeline_reports_the_missing_covariate_table():
    environment = os.environ.copy()
    environment.update({"COLUMNS": "80", "PYTHONDONTWRITEBYTECODE": "1"})
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "postgwas",
            "pipeline",
            "--modules",
            "magmacovar",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert completed.returncode == 2
    normalized = " ".join(completed.stdout.split())
    assert (
        "Required argument not provided: --covariates. Provide --covariates "
        "VALUE or set modules.magmacovar.input.covariates_file in the run "
        "configuration."
    ) in normalized
    required = {
        option.dest: option.config_path
        for option in REGISTRY.get("magmacovar").required_options
    }
    assert required["covariates"] == (
        "modules.magmacovar.input.covariates_file"
    )


def test_direct_magma_help_retains_direct_input_files():
    help_text = _help(("magma", "--help"), 120)
    normalized_help = " ".join(help_text.split())
    assert "--snp-location-file" in help_text
    assert "--p-value-file" in help_text
    assert "--resolve-variants-to-reference" in help_text
    assert "--duplicate-policy {err,lowest_p,remove}" in normalized_help
    assert "Default: lowest_p" in normalized_help
    assert "Available options: err, lowest_p, remove" in normalized_help
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


def test_aligned_help_wraps_marked_coloured_reference_blocks(monkeypatch):
    monkeypatch.setenv("COLUMNS", "50")
    parser = argparse.ArgumentParser(
        prog="postgwas demo",
        description=format_cli_help_block(
            "[bold cyan]Reference[/bold cyan]\n"
            "This deliberately long reference sentence must wrap to the active "
            "terminal width without exposing Rich markup."
        ),
        formatter_class=AlignedRichHelpFormatter,
    )

    help_text = parser.format_help()

    assert "[bold cyan]" not in help_text
    assert "Reference" in help_text
    assert all(len(line) <= 50 for line in help_text.splitlines())


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
        get_annot_ldblock_parser(), "ld_block_populations",
    ) == argparse.SUPPRESS
    assert action_default(
        get_annot_ldblock_parser(), "ld_region_dir",
    ) == argparse.SUPPRESS
    assert [population.value for population in annotation.populations] == [
        "EUR", "AFR", "EAS",
    ]

    clumping = configuration.modules.ld_clumping
    clumping_parser = get_ld_clump_parser()
    for destination in (
        "clumping_methods", "lead_p", "candidate_p", "r2_clump", "r2_lead",
        "ld_window_kb", "merge_dist", "missing_index_action",
        "missing_chromosome_action", "remove_mhc",
    ):
        assert action_default(clumping_parser, destination) == argparse.SUPPRESS
    assert action_default(
        get_ld_clumping_population_parser(), "population",
    ) == argparse.SUPPRESS
    assert action_default(
        get_pipeline_genome_build_parser(), "genome_build",
    ) == argparse.SUPPRESS
    shared_build_factory = "postgwas.cli.common:get_pipeline_genome_build_parser"
    assert shared_build_factory not in REGISTRY.get("annot_ldblock").parser_factories
    assert REGISTRY.get("annot_ldblock").genome_build_config_path is None
    annotation_required = {
        option.dest: option.config_path
        for option in REGISTRY.get("annot_ldblock").required_options
    }
    assert annotation_required["ld_region_dir"] == (
        "modules.ld_annotation.inputs.ld_region_dir"
    )
    for module in ("ld_clump", "mixer"):
        assert shared_build_factory in REGISTRY.get(module).parser_factories
    assert shared_build_factory not in REGISTRY.get("qc_summary").parser_factories
    assert REGISTRY.get("qc_summary").genome_build_config_path is None

    imputation = configuration.modules.imputation
    imputation_parser = get_common_imputation_parser()
    pred_ld = imputation.engines.pred_ld
    for destination in (
        "imputation_engine",
        "imputation_ld_reference",
        "imputation_r2_threshold",
        "imputation_minimum_maf",
        "ref",
        "corr_method",
        "resource_directory",
    ):
        assert action_default(imputation_parser, destination) == argparse.SUPPRESS
    assert imputation_parser.parse_args(
        ["--correlation-method", "spearman"],
    ).corr_method == "spearman"
    assert action_default(
        get_imputation_population_parser(), "population",
    ) == argparse.SUPPRESS
    normalized_imputation_help = " ".join(
        build_imputation_parser().format_help().split()
    )
    for expected in (
        "Default: pred_ld",
        "Default: %s" % pred_ld.minimum_r2,
        "Default: %s" % pred_ld.minimum_maf,
        "Default: %s" % pred_ld.mode,
        "Default: %s" % pred_ld.correlation_method,
        "Default: %s" % imputation.population.value,
        "Default: unset",
    ):
        assert expected in normalized_imputation_help

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
        assert action_default(manhattan_parser, destination) == argparse.SUPPRESS

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

    inferred = argparse.Namespace()
    pipeline_cli._resolve_pipeline_genome_build(
        inferred, ("ld_clump", "qc_summary"), configuration,
    )
    assert inferred.genome_build == "GRCh37"

    explicit = argparse.Namespace(genome_build="GRCh38")
    pipeline_cli._resolve_pipeline_genome_build(
        explicit, ("annot_ldblock", "ld_clump"), configuration,
    )
    assert explicit.genome_build == "GRCh38"


def test_ld_annotation_help_derives_build_from_vcf_and_colors_bed_warning():
    parser = build_ld_annotation_parser()
    destinations = {action.dest for action in parser._actions}
    assert "genome_build" not in destinations
    assert "genome_build" not in type(
        load_configuration().modules.ld_annotation
    ).model_fields

    action = next(
        action for action in parser._actions if action.dest == "ld_region_dir"
    )
    rendered = Text.from_markup(action.help or "")
    assert "Warning:" in rendered.plain
    assert "same genome build declared inside the input GWAS-VCF" in rendered.plain
    assert "<BUILD>_<POPULATION>_ldetect.bed.gz" in rendered.plain
    assert "GRCh37_EUR_ldetect.bed.gz" in rendered.plain
    assert "spelling and .bed.gz suffix must be exact" in rendered.plain
    assert "Required BED columns (tab-separated):" in rendered.plain
    assert "1. CHROM — exact input-VCF contig name" in rendered.plain
    assert "2. START — zero-based block start; included" in rendered.plain
    assert "3. END — zero-based block end; excluded" in rendered.plain
    assert "4. BLOCK_LABEL — nonempty value written to <POP>_LDblock" in rendered.plain
    assert "must end with _<START>_<END>" in rendered.plain
    assert "Additional trailing columns are allowed but ignored" in rendered.plain
    assert "Example row: 1<TAB>16103<TAB>2047216<TAB>EUR-1_16103_2047216" in (
        rendered.plain
    )
    assert "suffix must be exact.\n\nRequired BED columns" in rendered.plain
    assert "EUR-1_16103_2047216\n\nWarning:" in rendered.plain
    styles = [str(span.style) for span in rendered.spans]
    assert "bold cyan" in styles
    assert "bold bright_yellow" in styles
    assert "bright_yellow" in styles

    normalized = " ".join(_help(("annot_ldblock", "--help"), 120).split())
    assert "--genome-build" not in normalized
    assert "Genome coordinates:" not in normalized
    assert "Warning:" in normalized
    assert "<BUILD>_<POPULATION>_ldetect.bed.gz" in normalized
    assert "GRCh37_EUR_ldetect.bed.gz" in normalized
    assert "Required BED columns (tab-separated):" in normalized
    assert "Additional trailing columns are allowed but ignored" in normalized

    pipeline_help = " ".join(
        _help(("pipeline", "--modules", "annot_ldblock", "--help"), 120).split()
    )
    assert "--genome-build" not in pipeline_help
    assert "Genome coordinates:" not in pipeline_help
    assert "Warning:" in pipeline_help
    assert "<BUILD>_<POPULATION>_ldetect.bed.gz" in pipeline_help
    assert "GRCh37_EUR_ldetect.bed.gz" in pipeline_help
    assert "Required BED columns (tab-separated):" in pipeline_help
    assert "Additional trailing columns are allowed but ignored" in pipeline_help


@pytest.mark.parametrize(
    "arguments",
    (
        ("annot_ldblock", "--help"),
        ("pipeline", "--modules", "annot_ldblock", "--help"),
    ),
)
@pytest.mark.parametrize("width", (80, 120))
def test_ld_annotation_bed_contract_is_visually_separated(arguments, width):
    help_text = _help(arguments, width)
    assert re.search(
        r"suffix must be exact\.\n\s*\n\s+Required BED columns \(tab-separated\):",
        help_text,
    )
    assert re.search(
        r"EUR-1_16103_2047216\n\s*\n\s+Warning:",
        help_text,
    )


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
    assert "--samp-prev-warning-threshold" not in direct_help

    assert "Optional sample-prevalence override" in pipeline_help
    assert "PostGWAS uses the formatter-returned value" in pipeline_help
    assert "it does not recalculate it" in pipeline_help
    assert "--samp-prev-warning-threshold VALUE" in pipeline_help
    assert "Absolute difference above which" in pipeline_help
    assert "values are expressed as proportions" in pipeline_help
    assert "Default: 0.01" in pipeline_help

    for help_text in (direct_help, pipeline_help):
        assert help_text.count("Default: unset") >= 2


@pytest.mark.parametrize("width", (80, 120))
def test_ldsc_help_uses_magma_style_structure_and_alignment(width):
    direct_help = _help(("heritability", "--help"), width)
    pipeline_help = _help(
        ("pipeline", "--modules", "heritability", "--help"), width,
    )

    assert direct_help.startswith(
        "Usage: postgwas heritability --ldsc-input PATH --merge-alleles PATH"
    )
    for section in (
        "LDSC input:",
        "LDSC reference files:",
        "LDSC reference metadata:",
        "LDSC executables:",
        "Liability-scale settings:",
        "LDSC heritability settings:",
        "Munge-sumstats settings:",
    ):
        assert section in direct_help
    for label in (
        "Create the LDSC input first:",
        "Run observed-scale LDSC heritability:",
        "Run observed- and liability-scale LDSC heritability:",
        "Export reusable LDSC settings:",
        "Run with exported LDSC settings and explicit input resources:",
    ):
        assert label in direct_help
    assert "formatted/STUDY_ldsc_input.tsv" in direct_help
    assert "formatted/STUDY_ldsc_input.tsv.gz" not in direct_help

    for label in (
        "Run observed-scale LDSC heritability from a GWAS-VCF:",
        "Add liability-scale heritability using the GWAS-VCF case fraction:",
        "Override the GWAS-VCF case fraction for liability-scale heritability:",
    ):
        assert label in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help
    assert "--bcftools" not in pipeline_help
    assert "External software:" not in pipeline_help
    assert REGISTRY.get("heritability").pipeline_supplied_options == ()

    lines = direct_help.splitlines()
    population_line = next(
        line for line in lines if "--ldsc-population CODE" in line
    )
    population_column = population_line.index("Population represented")
    build_index = next(
        index for index, line in enumerate(lines)
        if "--ldsc-genome-build BUILD" in line
    )
    build_description = next(
        line for line in lines[build_index:build_index + 2]
        if "Declared genome build" in line
    )
    assert build_description.index("Declared genome build") == population_column


def test_ldsc_pipeline_examples_use_public_options_and_required_resources():
    option_strings = {
        option
        for action in build_ldsc_parser()._actions
        for option in action.option_strings
    }
    pipeline_only_options = {"--modules", "--vcf"}

    for _label, command, arguments in get_ldsc_pipeline_examples():
        assert command == "postgwas pipeline"
        example_options = {argument.split()[0] for argument in arguments}
        assert example_options <= option_strings | pipeline_only_options
        assert {
            "--modules",
            "--vcf",
            "--merge-alleles",
            "--ref-ld-chr",
            "--w-ld-chr",
            "--dataset-id",
            "--output-directory",
        } <= example_options


@pytest.mark.parametrize("width", (80, 160))
@pytest.mark.parametrize(
    "arguments",
    (
        ("heritability", "--help"),
        ("pipeline", "--modules", "heritability", "--help"),
    ),
)
def test_ldsc_help_exposes_only_not_m_5_50(arguments, width):
    help_text = _help(arguments, width)
    normalized = " ".join(help_text.split())

    assert "--use-M-5-50" not in help_text
    assert "--not-M-5-50" in help_text
    assert "this switch is not used by default" in normalized
    assert "Default: false" in normalized


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
        "❌  Pipeline failed: GCTA failed for [/analysis/pathway.set]"
    )
