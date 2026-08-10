"""Regression tests for the complete public command-line help surface."""

from __future__ import annotations

from io import StringIO
import os
import subprocess
import sys

import pytest
from rich.console import Console
from rich.text import Text

from postgwas.config import load_configuration
from postgwas.modules.enrichment.cli import get_geneset_parser
from postgwas.modules.caldera.cli import build_parser as build_caldera_parser
from postgwas.modules.flames.cli import build_parser as build_flames_parser
from postgwas.modules.gcta_gene.cli import get_gcta_gene_parser
from postgwas.modules.harmonisation.cli import get_harmonisation_parser
from postgwas.modules.magma.cli import build_parser as build_magma_parser
from postgwas.modules.magmacovar.cli import build_parser as build_magmacovar_parser
from postgwas.modules.mixer.cli import get_mixer_parser
from postgwas.modules.pops.cli import build_parser as build_pops_parser
from postgwas.modules.kpops.cli import build_parser as build_kpops_parser
from postgwas.cli.compute import get_compute_parser
from postgwas.core.ui import (
    format_cli_default,
    help_with_default,
    style_cli_defaults,
)
from postgwas.pipeline.registry import REGISTRY
from postgwas.pipeline import cli as pipeline_cli


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
    "sumstat_filter": ("sumstat_filter", "--help"),
    "validate": ("--validate", "--help"),
}

PIPELINE_TARGETS = tuple(
    name for name in REGISTRY.names() if REGISTRY.get(name).pipeline_enabled
)


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


@pytest.mark.parametrize("width", (80, 120))
def test_every_public_help_page_is_readable_and_has_examples(width):
    for command, arguments in PUBLIC_HELP_COMMANDS.items():
        help_text = _help(arguments, width)
        assert "postgwas" in help_text, command
        assert "example" in help_text.lower(), command
        assert "usage: postgwas-" not in help_text.lower(), command


def test_every_pipeline_target_builds_context_specific_help():
    for module in PIPELINE_TARGETS:
        help_text = _help(("pipeline", "--modules", module, "--help"), 120)
        assert "Usage: postgwas pipeline [options] --modules %s" % module in help_text
        assert "Workflow examples" in help_text


def test_pipeline_magma_help_hides_only_formatter_created_inputs():
    help_text = _help(("pipeline", "--modules", "magma", "--help"), 120)
    normalized_help = " ".join(help_text.split())

    assert "--resume, --no-resume" in help_text
    assert "--overwrite" in help_text
    assert "--resolve-variants-to-reference" in help_text
    assert "Configured default: false" in normalized_help
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


def test_pipeline_gcta_help_shows_analysis_options_and_specific_examples():
    help_text = _help(("pipeline", "--modules", "gcta_gene", "--help"), 120)

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
    assert "Run MAGMA gene analysis followed by PoPS prioritisation:" in help_text
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

    for option in (
        "--magma-association-prefix magma/results/STUDY_magma_35up_10down",
        "--kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv",
        "--kernel-matrix-prefix reference/kpops/kernel_linear",
        "--kpops-genome-build GRCh37",
    ):
        assert option in direct_help
    for option in (
        "--vcf study_GRCh37.vcf.gz",
        "--magma-ld-reference reference/1000G_EUR",
        "--gene-location-file reference/NCBI37.3.gene.loc",
        "--kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv",
        "--kernel-matrix-prefix reference/kpops/kernel_linear",
        "--kpops-genome-build GRCh37",
        "--dataset-id STUDY",
        "--output-directory results",
    ):
        assert option in pipeline_help
    assert "--run-config analysis_pipeline.yaml" not in pipeline_help


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
    assert "Configured default: false" in normalized_help


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


def test_default_label_and_value_use_distinct_central_styles():
    markup = format_cli_default(50)
    assert markup == (
        "[bold green]Default:[/bold green] [cyan]50[/cyan]"
    )
    rendered = Text.from_markup(markup)
    assert rendered.plain == "Default: 50"
    assert [str(span.style) for span in rendered.spans] == ["bold green", "cyan"]
    assert format_cli_default(False, label="Configured default") == (
        "[bold green]Configured default:[/bold green] [cyan]false[/cyan]"
    )
    assert help_with_default("Window size", "[configured]") == (
        "Window size ([bold green]Default:[/bold green] "
        "[cyan]\\[configured][/cyan])."
    )
    assert style_cli_defaults("Threshold (default: 0.2).") == (
        "Threshold ([bold green]Default:[/bold green] [cyan]0.2[/cyan])."
    )
    assert style_cli_defaults("Packaged default: EUR") == (
        "[bold green]Packaged default:[/bold green] [cyan]EUR[/cyan]"
    )
    assert style_cli_defaults("Configured YAML default: false") == (
        "[bold green]Configured YAML default:[/bold green] [cyan]false[/cyan]"
    )
    already_styled = "[bold green]Default:[/bold green] [cyan]0.2[/cyan]"
    assert style_cli_defaults(already_styled) == already_styled


def test_every_displayed_module_default_uses_shared_distinct_styles():
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
            assert "[bold green]" in help_text
            assert "[cyan]" in help_text
    assert displayed_defaults >= 20


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
