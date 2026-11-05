"""Documentation-command validation without executing scientific analyses."""

from __future__ import annotations

from pathlib import Path
import runpy
import shlex

import pytest


VALIDATOR = Path(__file__).parents[1] / "tools/docs/validate_wiki_cli.py"
MAGMA_REQUIRED_HELP = (
    "  --snp-location-file PATH  Required: SNP-coordinate input.\n"
    "  --p-value-file PATH\n"
    "                           Required: Variant P-value input.\n"
)


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
            "--modules gcta_gene --method fastbat_set --gmt '[cohort].gmt'",
            "--modules gcta_gene --method fastbat_set --gmt '[cohort].gmt'",
        ),
        (
            "--modules gcta_gene --method fastbat_set --gmt '[pipeline options].gmt'",
            "--modules gcta_gene --method fastbat_set --gmt '[pipeline options].gmt'",
        ),
        ("--modules flames '[pipeline options]'", "--modules flames"),
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


@pytest.mark.parametrize(("token", "expected"), [
    ("--analysis", "--analysis"),
    ("--analysis=gsa", "--analysis"),
    ("data --unknown file.tsv", None),
])
def test_only_actual_option_tokens_are_classified_as_options(validator, token, expected):
    assert validator["_clean_option"](token) == expected


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


@pytest.mark.parametrize("command,expected", [
    ("bash tools/setup/install_postgwas.sh --all-tools --update --name postgwas --jobs 2", ()),
    ("bash tools/setup/install_postgwas.sh --name 'postgwas' --work-directory '/scratch/my runs'", ()),
    ("docker run --name postgwas --rm postgwas:local postgwas --help", ("postgwas --help",)),
])
def test_environment_or_container_names_do_not_become_cli_invocations(
        validator, monkeypatch, tmp_path, command, expected):
    readme = tmp_path / "README.md"
    readme.write_text("```bash\n" + command + "\n```\n", encoding="utf-8")
    globals_ = validator["documented_commands"].__globals__
    monkeypatch.setitem(globals_, "ROOT_README", readme)
    monkeypatch.setitem(globals_, "_published_sources", lambda manifest: ())
    assert validator["documented_commands"](include_inline=False) == tuple(
        (readme, invocation) for invocation in expected)


def test_validator_uses_selected_help_and_separates_cache_entries(
    validator, monkeypatch,
):
    commands = [
        "postgwas pipeline --modules finemap --help",
        "postgwas pipeline --modules finemap --clumping-methods standard "
        "--finemap-method susie --vcf study.vcf.gz --help",
        "postgwas pipeline --modules finemap --clumping-methods standard "
        "--finemap-method susie --vcf another_study.vcf.gz --help",
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
        lambda manifest, **kwargs: tuple((VALIDATOR, command) for command in commands),
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
        lambda manifest, **kwargs: ((VALIDATOR, "postgwas pipeline --modules mixer --analysis wrong"),),
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


@pytest.mark.parametrize("usage", [
    "Usage: postgwas magma --snp-location-file PATH --p-value-file PATH [--gene-set-file PATH]\n\nOptions:",
    "usage: postgwas magma --snp-location-file PATH\n                      --p-value-file PATH [options]\n\nOptions:",
])
def test_required_usage_ignores_optional_and_wrapped_groups(validator, usage):
    help_text = usage + "\n" + MAGMA_REQUIRED_HELP + (
        "  --gene-set-file PATH      Required: When running gene-set analysis.\n"
    )
    assert validator["_usage_required_options"](help_text) == {"--snp-location-file", "--p-value-file"}


def test_incomplete_magma_recipe_is_rejected(validator, monkeypatch):
    globals_ = validator["validate_documented_commands"].__globals__
    monkeypatch.setitem(globals_, "_command_inventory", lambda executable: ({"magma"}, set()))
    monkeypatch.setitem(globals_, "documented_commands", lambda manifest, **kwargs: (
        (VALIDATOR, "postgwas magma --dataset-id STUDY"),
    ))
    monkeypatch.setitem(globals_, "_run_help", lambda executable, arguments:
        "Usage: postgwas magma --snp-location-file PATH --p-value-file PATH [options]\n\n"
        + MAGMA_REQUIRED_HELP +
        "  --dataset-id NAME\n")
    with pytest.raises(validator["CliDocumentationError"], match="incomplete CLI recipe") as error:
        validator["validate_documented_commands"]()
    assert "--p-value-file, --snp-location-file" in str(error.value)


@pytest.mark.parametrize("option", ["--p-value-file values.tsv", "--p-value-file=values.tsv"])
def test_complete_recipe_accepts_both_value_syntaxes(validator, option):
    tokens = shlex.split("postgwas magma --snp-location-file snps.tsv " + option)
    assert not validator["_missing_recipe_options"](tokens,
        "Usage: postgwas magma --snp-location-file PATH --p-value-file PATH [options]\n\n"
        + MAGMA_REQUIRED_HELP)


def test_recipe_bare_required_flag_is_not_a_value(validator):
    assert validator["_missing_recipe_options"](
        ["postgwas", "magma", "--p-value-file"],
        "Usage: postgwas magma --p-value-file PATH\n\n" + MAGMA_REQUIRED_HELP,
    ) == {"--p-value-file"}


def test_placeholder_yaml_and_help_are_not_claimed_to_be_complete_recipes(validator):
    for args in (["--run-config", "resources.yaml"], ["--help"]):
        assert not validator["_missing_recipe_options"](
            ["postgwas", "magma", *args], "Usage: postgwas magma --p-value-file PATH\n\n"
            + MAGMA_REQUIRED_HELP)


def test_defaulted_direct_susie_selector_is_not_an_invented_requirement(validator):
    help_text = (
        "Usage: postgwas finemap --finemap-method {susie,finemap} [options]\n\n"
        "  --finemap-method METHOD  Fine-mapping method. (Default: susie)\n"
        "                           Available options: susie, finemap.\n"
        "  --susie-input-file PATH  Required: For the SuSiE method.\n"
        "  --locus-file PATH        Required: Locus table.\n"
        "  --finemap-ld-reference PREFIX  Required: PLINK reference.\n"
    )
    tokens = shlex.split(
        "postgwas finemap --susie-input-file 'study [meta].tsv' --locus-file loci.tsv "
        "--finemap-ld-reference reference --dataset-id STUDY --output-directory results"
    )
    assert validator["_usage_required_options"](help_text) == set()
    assert validator["_missing_recipe_options"](tokens, help_text) == set()
    assert validator["_choice_options"](help_text)["--finemap-method"] == {"susie", "finemap"}


def test_required_option_alias_and_wrapped_label_are_preserved(validator):
    help_text = (
        "Usage: postgwas magma --p-value-file PATH [options]\n\n"
        "  --p-value-file PATH, --pval-file PATH\n"
        "                            Required: Variant P-values.\n"
        "  --reference PATH          This mentions Required: in ordinary prose.\n"
    )
    assert validator["_usage_required_options"](help_text) == {"--p-value-file"}
    assert not validator["_missing_recipe_options"](
        shlex.split("postgwas magma --pval-file values.tsv"), help_text,
    )


@pytest.mark.parametrize("arguments", [
    ["[options]"], ["[overrides]"], ["[arguments]"], ["[args]"],
    ["[pipeline", "options]"], ["[pipeline options]"],
])
def test_executable_synopsis_abbreviations_are_rejected(validator, arguments):
    with pytest.raises(validator["CliDocumentationError"], match="must not abbreviate"):
        validator["_missing_recipe_options"](["postgwas", "magma", *arguments], "")


@pytest.mark.parametrize("filename", [
    "study[meta].tsv", "data/[cohort]/study.tsv", "data [meta]/study.tsv",
    "data --help file.tsv", "data --run-config file.tsv", "data --unknown file.tsv",
])
def test_harvested_quoted_filenames_preserve_tokens_and_required_checks(
    validator, monkeypatch, tmp_path, filename,
):
    tokens = ["postgwas", "magma", "--snp-location-file", filename]
    readme = tmp_path / "README.md"
    readme.write_text("```bash\n" + shlex.join(tokens) + "\n```\n", encoding="utf-8")
    globals_ = validator["documented_commands"].__globals__
    monkeypatch.setitem(globals_, "ROOT_README", readme)
    monkeypatch.setitem(globals_, "_published_sources", lambda manifest: ())
    commands = validator["documented_commands"](include_inline=False)
    assert len(commands) == 1
    harvested = shlex.split(commands[0][1])
    assert harvested == tokens
    help_text = (
        "Usage: postgwas magma --snp-location-file PATH --p-value-file PATH [options]\n\n"
        + MAGMA_REQUIRED_HELP
    )
    assert validator["_missing_recipe_options"](harvested, help_text) == {"--p-value-file"}
    monkeypatch.setitem(globals_, "REPOSITORY_ROOT", tmp_path)
    monkeypatch.setitem(globals_, "_command_inventory", lambda executable: ({"magma"}, set()))
    monkeypatch.setitem(globals_, "_run_help", lambda executable, arguments: help_text)
    with pytest.raises(validator["CliDocumentationError"], match="incomplete CLI recipe"):
        validator["validate_documented_commands"]()


def test_pipeline_recipe_checks_upstream_references_not_generated_inputs(validator):
    required = validator["_pipeline_required_options"](
        ["postgwas", "pipeline", "--modules", "magma"])
    assert {"--vcf", "--dataset-id", "--output-directory", "--magma-ld-reference", "--gene-location-file"} <= required
    assert not {"--snp-location-file", "--p-value-file", "--format"} & required


def test_fenced_and_inline_commands_have_distinct_completeness_scope(validator, monkeypatch, tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text("`postgwas magma`\n```bash\npostgwas magma --help\n```\n", encoding="utf-8")
    globals_ = validator["documented_commands"].__globals__
    monkeypatch.setitem(globals_, "ROOT_README", readme)
    monkeypatch.setitem(globals_, "_published_sources", lambda manifest: ())
    assert validator["documented_commands"](include_inline=False) == ((readme, "postgwas magma --help"),)
    assert len(validator["documented_commands"]()) == 2
