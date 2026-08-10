"""Tests for deterministic publication of canonical docs to GitHub Wiki."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import yaml


REPOSITORY_ROOT = Path(__file__).parents[1]
ROOT_README = REPOSITORY_ROOT / "README.md"
BUILDER = REPOSITORY_ROOT / "tools" / "docs" / "build_wiki.py"
CLI_VALIDATOR = REPOSITORY_ROOT / "tools" / "docs" / "validate_wiki_cli.py"
MODULE_TEMPLATE = REPOSITORY_ROOT / "docs" / "templates" / "module-page.md"
MANIFEST = REPOSITORY_ROOT / "docs" / "wiki.yml"
MODULE_DIRECTORY = REPOSITORY_ROOT / "docs" / "wiki" / "modules"
HARMONISATION_ORDER_SOURCE = (
    REPOSITORY_ROOT / "docs" / "wiki" / "harmonisation" / "processing-order.md"
)
HARMONISATION_ORDER_ASSET = HARMONISATION_ORDER_SOURCE.with_name(
    "harmonisation-processing-order.svg"
)
CONFIGURATION_DEFAULTS_SOURCE = (
    REPOSITORY_ROOT / "docs" / "wiki" / "reference" / "configuration-defaults.md"
)
HARMONISATION_CONFIGURATION_SOURCE = (
    REPOSITORY_ROOT / "docs" / "modules" / "harmonisation" / "configuration.md"
)
HARMONISATION_SAMPLE_SHEET_SOURCE = (
    REPOSITORY_ROOT / "docs" / "wiki" / "harmonisation" / "sample-sheet.md"
)
WIKI_WORKFLOW = REPOSITORY_ROOT / ".github" / "workflows" / "wiki-docs.yml"
MODULE_HEADINGS = [
    "## Purpose",
    "## What the analysis does",
    "## When to use it",
    "## Input requirements",
    "## Command",
    "## Minimal example",
    "## Full example",
    "## Parameters",
    "## Processing steps",
    "## Outputs",
    "## QC and logs",
    "## Interpretation",
    "## Common problems",
    "## Limitations",
    "## Scientific references",
]


def _wiki_pages():
    manifest = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    return manifest["wiki"]["pages"]


def test_module_template_preserves_the_researcher_user_journey():
    template = MODULE_TEMPLATE.read_text(encoding="utf-8")
    positions = [template.index(heading) for heading in MODULE_HEADINGS]
    assert positions == sorted(positions)
    assert "What does this do? → What do I need? → How do I run it? → What do I get?" in template


def test_every_analysis_module_follows_the_standard_page_structure():
    for path in sorted(MODULE_DIRECTORY.glob("*.md")):
        text = path.read_text(encoding="utf-8")
        actual = [line for line in text.splitlines() if line.startswith("## ")]
        assert actual == MODULE_HEADINGS, path


def test_wiki_sources_validate_without_writing(tmp_path):
    output = tmp_path / "wiki"
    completed = subprocess.run(
        [sys.executable, str(BUILDER), "--check", "--output", str(output)],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert f"Wiki validated: {len(_wiki_pages())} pages" in completed.stdout
    assert not output.exists()


def test_harmonisation_wiki_preserves_the_runtime_processing_order():
    text = HARMONISATION_ORDER_SOURCE.read_text(encoding="utf-8")
    dataset_steps = [
        "**Dataset step 1 — Check configuration.**",
        "**Dataset step 2 — Check the header.**",
        "**Dataset step 3 — Read and clean the study.**",
        "**Dataset step 4 — Check the parsed values.**",
        "**Dataset step 5 — Determine the genome build.**",
        "**Dataset step 6 — Determine strand consensus.**",
        "**Dataset step 7 — Determine statistic types.**",
        "**Dataset step 8 — Split by chromosome.**",
    ]
    chromosome_steps = [
        "**Chromosome step 1 — Load the partition.**",
        "**Chromosome step 2 — Find the reference files.**",
        "**Chromosome step 3 — Put effects on the BETA scale.**",
        "**Chromosome step 4 — Orient alleles and obtain EAF.**",
        "**Chromosome step 5 — Calculate sample size.**",
        "**Chromosome step 6 — Recover BETA or SE from Z when needed.**",
        "**Chromosome step 7 — Harmonise P values.**",
        "**Chromosome step 8 — Derive a still-missing SE.**",
        "**Chromosome step 9 — Produce the final Z score.**",
        "**Chromosome step 10 — Check statistical agreement.**",
        "**Chromosome step 11 — Obtain INFO.**",
        "**Chromosome step 12 — Complete variant identifiers.**",
        "**Chromosome step 13 — Check required output fields.**",
        "**Chromosome step 14 — Export adapter input.**",
        "**Chromosome step 15 — Create the first VCF.**",
        "**Chromosome step 16 — Normalize, annotate, and lift over.**",
    ]
    post_merge_steps = [
        "**Post-merge step 1 — Merge chromosome VCFs.**",
        "**Post-merge step 2 — Compare population frequencies.**",
        "**Post-merge step 3 — Combine side files and clean up.**",
        "**Post-merge step 4 — Assess the raw merged VCF.**",
        "**Post-merge step 5 — Finish the dataset.**",
    ]

    for sequence in (dataset_steps, chromosome_steps, post_merge_steps):
        positions = [text.index(heading) for heading in sequence]
        assert positions == sorted(positions)


def test_harmonisation_cross_page_step_and_info_contracts_match_runtime():
    configuration = HARMONISATION_CONFIGURATION_SOURCE.read_text(encoding="utf-8")
    sample_sheet = HARMONISATION_SAMPLE_SHEET_SOURCE.read_text(encoding="utf-8")

    assert "At chromosome step 04, the supplied raw frequency reference" in configuration
    assert "At\nchromosome step 03, an odds ratio is normalized" in configuration
    assert "reciprocates an odds ratio before its\nlater log conversion" not in configuration
    assert "Internal INFO takes priority when both are listed" in sample_sheet
    assert "Exactly one source is required" in sample_sheet


def test_readme_uses_current_repository_and_version_two_interfaces():
    text = ROOT_README.read_text(encoding="utf-8")

    assert "## Documentation" in text
    assert "GitHub Wiki access is not enabled" in text
    assert "JIBINJOHNV/postgwas_v2/wiki" not in text
    assert "JIBINJOHNV/postgwas.git" not in text
    assert "jibinjv/postgwas:1.3" not in text
    assert "sumstat_file" not in text
    assert "--sample-sheet studies.csv" in text


def test_readme_indexes_every_manifest_page_in_order():
    text = ROOT_README.read_text(encoding="utf-8")
    positions = []
    for page in _wiki_pages():
        link = f']({page["source"]})'
        assert link in text, page
        positions.append(text.index(link))

    assert positions == sorted(positions)


def test_documentation_workflow_validates_without_wiki_publication():
    workflow = WIKI_WORKFLOW.read_text(encoding="utf-8")

    assert workflow.count('- "README.md"') == 2
    assert "name: Validate user documentation" in workflow
    assert "python -m pytest -p no:cacheprovider tests/test_wiki_docs.py -q" in workflow
    assert "python tools/docs/validate_wiki_cli.py" in workflow
    assert "publish:" not in workflow
    assert ".wiki.git" not in workflow


def test_harmonisation_wiki_visualises_execution_scope_and_variant_outcomes():
    text = HARMONISATION_ORDER_SOURCE.read_text(encoding="utf-8")

    assert "```mermaid" not in text
    assert "(harmonisation-processing-order.svg)" in text
    root = ET.fromstring(
        HARMONISATION_ORDER_ASSET.read_text(encoding="utf-8")
    )
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    assert root.find("svg:title", namespace) is not None
    assert root.find("svg:desc", namespace) is not None
    labels = " ".join(root.itertext())
    assert "Once for the complete study" in labels
    assert "Repeated for every chromosome" in labels
    assert "Swapped or complement-swapped" in labels
    assert "Ambiguous or unmatched" in labels
    assert "Retry failed chromosome only" in labels
    assert "Optional input-to-VCF validation" in labels


def test_published_wiki_uses_supported_command_names_and_config_actions():
    completed = subprocess.run(
        [sys.executable, str(CLI_VALIDATOR), "--commands-only"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "documented commands (command names)" in completed.stdout


def test_wiki_build_creates_navigation_and_rewrites_source_links(tmp_path):
    output = tmp_path / "wiki"
    subprocess.run(
        [sys.executable, str(BUILDER), "--output", str(output)],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    expected_pages = {f'{page["slug"]}.md' for page in _wiki_pages()}
    assert {path.name for path in output.iterdir()} == expected_pages | {
        ".postgwas-wiki-files",
        "_Footer.md",
        "_Sidebar.md",
        HARMONISATION_ORDER_ASSET.name,
    }
    assert (output / HARMONISATION_ORDER_ASSET.name).read_bytes() == (
        HARMONISATION_ORDER_ASSET.read_bytes()
    )
    home = (output / "Home.md").read_text(encoding="utf-8")
    assert "Generated from docs/wiki/home.md" in home
    assert "[Configuration](Configuration)" in home
    assert "core/configuration.md" not in home

    sidebar = (output / "_Sidebar.md").read_text(encoding="utf-8")
    assert sidebar.startswith("# PostGWAS User Guide\n")
    assert "- [Home](Home)" in sidebar
    assert "## Core Concepts" in sidebar
    assert "- [Logging and Reproducibility](Logging-and-Reproducibility)" in sidebar

    defaults_source = CONFIGURATION_DEFAULTS_SOURCE.read_text(encoding="utf-8")
    defaults_page = (output / "Configuration-Defaults.md").read_text(encoding="utf-8")
    assert "<!-- GENERATED: PACKAGED CONFIGURATION DEFAULTS -->" in defaults_source
    assert "<!-- GENERATED: PACKAGED CONFIGURATION DEFAULTS -->" not in defaults_page
    assert "maf_min: 0.01" in defaults_page
    assert "modules:\n  allele_orientation:" in defaults_page
    assert "  harmonisation:\n    enabled: false" in defaults_page
    assert "scientific_policy:" not in defaults_page


def test_wiki_build_is_deterministic(tmp_path):
    output = tmp_path / "wiki"
    command = [sys.executable, str(BUILDER), "--output", str(output)]
    subprocess.run(command, cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True)
    first = {
        path.name: path.read_bytes()
        for path in output.iterdir()
    }
    subprocess.run(command, cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True)
    second = {
        path.name: path.read_bytes()
        for path in output.iterdir()
    }

    assert second == first


def test_wiki_rebuild_preserves_unmanaged_files(tmp_path):
    output = tmp_path / "wiki"
    command = [sys.executable, str(BUILDER), "--output", str(output)]
    subprocess.run(command, cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True)
    unrelated = output / "user-owned.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    subprocess.run(command, cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True)

    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_wiki_manifest_rejects_duplicate_slugs(tmp_path):
    manifest = tmp_path / "wiki.yml"
    manifest.write_text(
        """\
wiki:
  title: Invalid Wiki
  footer_source: docs/wiki/footer.md
  pages:
    - title: Home
      slug: Home
      source: docs/wiki/home.md
    - title: Duplicate
      slug: home
      source: docs/architecture.md
""",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [sys.executable, str(BUILDER), "--manifest", str(manifest), "--check"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    assert "Duplicate Wiki slug: home" in completed.stderr
