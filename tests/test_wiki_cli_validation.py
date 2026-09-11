"""Documentation-command validation without executing scientific analyses."""

from __future__ import annotations

from pathlib import Path
import runpy
import shlex

import pytest


VALIDATOR = Path(__file__).parents[1] / "tools/docs/validate_wiki_cli.py"


@pytest.fixture
def validator():
    return runpy.run_path(str(VALIDATOR))


@pytest.mark.parametrize(
    ("arguments", "expected"),
    [
        ("", ""),
        ("--modules", ""),
        ("--modules finemap", "--modules finemap"),
        ("--modules flames [pipeline options]", "--modules flames"),
        (
            "--modules finemap --clumping-methods standard --finemap-method susie "
            "--vcf study.vcf.gz --run-config missing.yaml --threads 2",
            "--modules finemap --clumping-methods standard --finemap-method susie",
        ),
        (
            "--modules=finemap --clumping-methods=standard cojo-slct "
            "--finemap-method=finemap",
            "--modules finemap --clumping-methods standard cojo-slct "
            "--finemap-method finemap",
        ),
        (
            "--modules single_cell --tools scdrs ldsc_celltype --vcf study.vcf.gz",
            "--modules single_cell --tools scdrs ldsc_celltype",
        ),
        ("--modules mixer --analysis gsa", "--modules mixer --analysis gsa"),
        (
            "--modules gcta_cojo --cojo-mode cond",
            "--modules gcta_cojo --cojo-mode cond",
        ),
        (
            "--modules gcta_gene --method fastbat_set --gmt pathways.gmt",
            "--modules gcta_gene --method fastbat_set --gmt pathways.gmt",
        ),
        (
            "--modules gcta_gene --method fastbat_set --fastbat-set-list sets.txt",
            "--modules gcta_gene --method fastbat_set --fastbat-set-list sets.txt",
        ),
        (
            "--modules magma --apply-filter --apply-imputation --apply-manhattan "
            "--heritability",
            "--modules magma --apply-filter --apply-imputation --apply-manhattan "
            "--heritability",
        ),
        ("--apply-filter", "--apply-filter"),
    ],
)
def test_pipeline_help_preserves_documented_workflow_context(
    validator, arguments, expected,
):
    tokens = shlex.split("postgwas pipeline " + arguments)
    assert validator["_pipeline_help_arguments"](tokens) == [
        "pipeline", *shlex.split(expected), "--help",
    ]


@pytest.mark.parametrize(
    "selection",
    [
        "--modules unknown_target",
        "--modules finemap --clumping-methods wrong --finemap-method wrong",
        "--modules single_cell --tools wrong",
        "--modules mixer --analysis wrong",
        "--modules gcta_gene --method wrong",
        "--modules gcta_cojo --cojo-mode wrong",
    ],
)
def test_invalid_selectors_are_not_discarded_or_replaced(validator, selection):
    tokens = shlex.split("postgwas pipeline " + selection)
    assert validator["_pipeline_help_arguments"](tokens) == [
        "pipeline", *shlex.split(selection), "--help",
    ]


@pytest.mark.parametrize("option", ["--tools", "--method=", "--finemap-method"])
def test_missing_selector_value_is_not_silently_discarded(validator, option):
    tokens = shlex.split("postgwas pipeline --modules finemap " + option)
    with pytest.raises(validator["CliDocumentationError"], match="has no value"):
        validator["_pipeline_help_arguments"](tokens)


@pytest.mark.parametrize(
    "help_text",
    [
        "  --analysis {univariate,gsa,all}  Select analysis.\n",
        "  --analysis ANALYSIS  Choices: {univariate,gsa,all}.\n",
        "  --analysis ANALYSIS  Select analysis. Available options: univariate, gsa, all.\n",
        "  --analysis ANALYSIS  Select analysis. Available options:\n"
        "                       univariate, gsa, all.\n",
    ],
)
def test_choice_reader_accepts_argparse_and_shared_help_styles(validator, help_text):
    assert validator["_choice_options"](help_text) == {
        "--analysis": {"univariate", "gsa", "all"},
    }


def test_choice_reader_does_not_treat_description_references_as_option_aliases(validator):
    help_text = (
        "  --duplicate-policy POLICY  Resolve duplicates before\n"
        "                            --merge-alleles. Available options: error, exclude_all.\n"
        "  --merge-alleles PATH       Read the SNP-and-allele table.\n"
        "  --method METHOD           Select the method; --gmt supplies pathways.\n"
        "                            --fastbat-set-list is another source.\n"
        "                            Available options: fastbat_gene, fastbat_set.\n"
        "  --gmt PATH                Read the pathway table.\n"
        "  --reference PATH          A reference template such as chr{chromosome}.\n"
    )
    assert validator["_choice_options"](help_text) == {
        "--duplicate-policy": {"error", "exclude_all"},
        "--method": {"fastbat_gene", "fastbat_set"},
    }


def test_choice_reader_preserves_real_option_aliases(validator):
    assert validator["_choice_options"](
        "  --analysis MODE, --mode MODE  Available options: univariate, gsa.\n"
    ) == {
        "--analysis": {"univariate", "gsa"},
        "--mode": {"univariate", "gsa"},
    }


@pytest.mark.parametrize("selection", ["--analysis invalid", "--analysis=invalid"])
def test_invalid_choice_is_rejected_in_both_option_syntaxes(validator, selection):
    errors = validator["_documented_choice_errors"](
        shlex.split("postgwas pipeline --modules mixer " + selection),
        {"--analysis": {"univariate", "gsa", "all"}},
    )
    assert len(errors) == 1
    assert "--analysis value(s) invalid" in errors[0]


@pytest.mark.parametrize("selection", ["--analysis gsa", "--analysis=gsa"])
def test_valid_choice_is_accepted_in_both_option_syntaxes(validator, selection):
    assert validator["_documented_choice_errors"](
        shlex.split("postgwas pipeline --modules mixer " + selection),
        {"--analysis": {"univariate", "gsa", "all"}},
    ) == []


def test_current_installation_commands_are_not_harvested_as_analyses(
    validator, monkeypatch, tmp_path,
):
    readme = tmp_path / "README.md"
    readme.write_text(
        "```bash\n"
        "bash tools/setup/install_postgwas.sh --all-tools --name postgwas\n"
        "conda activate postgwas\n"
        "postgwas --help\n"
        "```\n",
        encoding="utf-8",
    )
    globals_ = validator["documented_commands"].__globals__
    monkeypatch.setitem(globals_, "ROOT_README", readme)
    monkeypatch.setitem(globals_, "_published_sources", lambda manifest: ())
    assert validator["documented_commands"]() == ((readme, "postgwas --help"),)


def test_validator_uses_selected_help_and_separates_cache_entries(
    validator, monkeypatch,
):
    commands = [
        "postgwas pipeline --modules finemap --help",
        "postgwas pipeline --modules finemap --clumping-methods standard "
        "--finemap-method susie --vcf study.vcf.gz",
        "postgwas pipeline --modules finemap --clumping-methods standard "
        "--finemap-method susie --vcf another_study.vcf.gz",
    ]
    calls = []

    def run_help(executable, arguments):
        calls.append(arguments)
        text = "  --modules MODULE\n  --clumping-methods METHOD\n  --finemap-method METHOD\n"
        if "--finemap-method" in arguments:
            text += "  --vcf PATH\n"
        return text

    globals_ = validator["validate_documented_commands"].__globals__
    monkeypatch.setitem(globals_, "_command_inventory", lambda executable: ({"pipeline"}, set()))
    monkeypatch.setitem(
        globals_, "documented_commands",
        lambda manifest: tuple((VALIDATOR, command) for command in commands),
    )
    monkeypatch.setitem(globals_, "_run_help", run_help)
    assert validator["validate_documented_commands"]() == 3
    assert calls == [
        ["pipeline", "--modules", "finemap", "--help"],
        ["pipeline", "--modules", "finemap", "--clumping-methods", "standard",
         "--finemap-method", "susie", "--help"],
    ]


def test_validator_reports_invalid_current_style_choice(validator, monkeypatch):
    globals_ = validator["validate_documented_commands"].__globals__
    monkeypatch.setitem(globals_, "_command_inventory", lambda executable: ({"pipeline"}, set()))
    monkeypatch.setitem(
        globals_, "documented_commands",
        lambda manifest: ((VALIDATOR, "postgwas pipeline --modules mixer --analysis wrong"),),
    )
    monkeypatch.setitem(
        globals_, "_run_help",
        lambda executable, arguments: (
            "  --modules MODULE\n"
            "  --analysis ANALYSIS  Available options: univariate, gsa, all.\n"
        ),
    )
    with pytest.raises(validator["CliDocumentationError"], match="unsupported choice"):
        validator["validate_documented_commands"]()
