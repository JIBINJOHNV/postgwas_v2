"""Tests for deterministic publication of canonical docs to GitHub Wiki."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import yaml

from postgwas.config import load_configuration
from postgwas.config.models.modules.single_cell import single_cell_pipeline_dependencies
from postgwas.modules.harmonisation.sample_sheet import HarmonisationSampleSheetRow
from postgwas.pipeline.planner import (
    build_pipeline_plan,
    resolve_pipeline_dependency_overrides,
)
from postgwas.pipeline.registry import REGISTRY


REPOSITORY_ROOT = Path(__file__).parents[1]
ROOT_README = REPOSITORY_ROOT / "README.md"
BUILDER = REPOSITORY_ROOT / "tools" / "docs" / "build_wiki.py"
CLI_VALIDATOR = REPOSITORY_ROOT / "tools" / "docs" / "validate_wiki_cli.py"
POLICY_SYNCHRONIZER = (
    REPOSITORY_ROOT / "tools" / "docs" / "update_harmonisation_policies.py"
)
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
HARMONISATION_POLICY_SOURCE = (
    REPOSITORY_ROOT / "docs" / "wiki" / "harmonisation" / "policies.md"
)
DOCUMENTATION_HOME = REPOSITORY_ROOT / "docs" / "wiki" / "home.md"
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
        "**Dataset step 6 — Determine strand policy and consensus.**",
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

    assert "## Documentation and support" in text
    assert "](docs/wiki/home.md)" in text
    assert "JIBINJOHNV/postgwas_v2/wiki" not in text
    assert "JIBINJOHNV/postgwas.git" not in text
    assert "jibinjv/postgwas:1.3" not in text
    assert "sumstat_file" not in text
    assert "--sample-sheet studies.csv" in text


def test_harmonisation_policy_guide_preserves_the_complete_runtime_order():
    section = HARMONISATION_POLICY_SOURCE.read_text(encoding="utf-8")

    stage_headings = [
        "## Stage 1 — preflight before the large input is read",
        "## Stage 2 — eight complete-dataset steps",
        "## Stage 3 — sixteen steps repeated for every chromosome",
        "## Stage 4 — chromosome retry and row reconciliation",
        "## Stage 5 — five post-merge steps",
        "## Stage 6 — optional input-to-VCF concordance",
        "## Stage 7 — multi-dataset summary and final states",
    ]
    positions = [section.index(heading) for heading in stage_headings]
    assert positions == sorted(positions)

    def numbered_rows(first_heading: str, second_heading: str) -> list[int]:
        table = section[
            section.index(first_heading) : section.index(second_heading)
        ]
        return [
            int(value)
            for value in re.findall(r"^\| (\d+) \|", table, flags=re.MULTILINE)
        ]

    assert numbered_rows(stage_headings[1], stage_headings[2]) == list(range(1, 9))
    assert numbered_rows(stage_headings[2], stage_headings[3]) == list(range(1, 17))
    assert numbered_rows(stage_headings[4], stage_headings[5]) == list(range(1, 6))


def test_harmonisation_policy_guide_preserves_every_sample_sheet_field():
    text = HARMONISATION_POLICY_SOURCE.read_text(encoding="utf-8")
    start = text.index("<summary>Complete version-2 sample-sheet field contract</summary>")
    contract = text[start : text.index("</details>", start)]

    for field in HarmonisationSampleSheetRow.model_fields:
        assert f"`{field}`" in contract, field


def test_harmonisation_policy_reference_matches_canonical_yaml():
    completed = subprocess.run(
        [sys.executable, str(POLICY_SYNCHRONIZER), "--check"],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Harmonisation policy reference is synchronized." in completed.stdout


def test_harmonisation_se_tail_documentation_matches_missing_cell_recovery():
    source = (
        REPOSITORY_ROOT / "src" / "postgwas" / "config" / "defaults"
        / "modules" / "harmonisation.yaml"
    )
    policy = yaml.safe_load(source.read_text(encoding="utf-8"))["policies"]["pvalue"]["se_tail"]
    text = " ".join(policy["help"].split())
    assert "chromosome step 08" in text
    assert "unresolved null cells in an existing SE column" in text
    assert "pvalue.derive_partial_missing_se is true" in text
    assert "study-supplied and Z-derived SE values are preserved" in text
    assert "pvalue.zero_missing_se separately" in text
    assert "only when the study supplies no standard-error column" not in text


def test_documentation_home_indexes_every_manifest_page_in_order():
    text = DOCUMENTATION_HOME.read_text(encoding="utf-8")
    documentation = text[text.index("## Complete documentation index") :]
    positions = []
    for page in _wiki_pages():
        source = REPOSITORY_ROOT / page["source"]
        # Home is the full documentation index; paths are relative to that page.
        link = f']({Path(os.path.relpath(source, DOCUMENTATION_HOME.parent)).as_posix()})'
        assert link in documentation, page
        positions.append(documentation.index(link))

    assert positions == sorted(positions)


def test_readme_lists_public_commands_and_supported_execution_modes():
    text = ROOT_README.read_text(encoding="utf-8")

    for command in REGISTRY.commands():
        assert (
            f"`{command}`" in text or f"`postgwas {command}`" in text
        ), command

    for name in REGISTRY.names():
        specification = REGISTRY.get(name)
        if specification.pipeline_enabled and specification.runner:
            assert f"`{name}`" in text, name

    assert "Direct mode" in text
    assert "Pipeline mode" in text
    for command in ("harmonisation", "pathway_enrichment"):
        catalogue_row = next(
            line for line in text.splitlines()
            if line.startswith("|") and f"`{command}`" in line
        )
        assert "standalone" in catalogue_row.lower(), catalogue_row
    assert "`qc_summary`" in text
    assert "allele_orientation" not in text


def test_readme_prioritises_installation_and_routes_to_detailed_guides():
    text = ROOT_README.read_text(encoding="utf-8")
    headings = [
        "## Installation",
        "## Available analyses",
        "## Choose your starting point",
        "## Harmonise your summary statistics",
        "## Run downstream analyses",
        "## Reference data and configuration",
        "## Results, logging and reproducibility",
        "## Documentation and support",
        "## Development and licensing",
    ]
    positions = [text.index(heading) for heading in headings]
    assert positions == sorted(positions)
    assert "bash tools/setup/install_postgwas.sh --all-tools" in text
    assert "](docs/wiki/harmonisation/policies.md)" in text
    assert "](docs/wiki/home.md)" in text
    assert "<!-- BEGIN GENERATED HARMONISATION POLICY REFERENCE -->" not in text
    guide = HARMONISATION_POLICY_SOURCE.read_text(encoding="utf-8")
    assert guide.count("<!-- BEGIN GENERATED HARMONISATION POLICY REFERENCE -->") == 1
    assert guide.count("<!-- END GENERATED HARMONISATION POLICY REFERENCE -->") == 1


def test_readme_connects_preparation_execution_and_result_review():
    text = ROOT_README.read_text(encoding="utf-8")
    headings = (
        "## Workflow at a glance",
        "### Before your first analysis",
        "### Quick start: from raw statistics to results",
        "### Generate and review the sample sheet",
        "### Complete the study information",
        "### Pipeline targets and execution order",
        "### Optional pipeline stages",
        "### Check reference compatibility",
        "### Check that the run succeeded",
        "### Moving from the previous PostGWAS version",
    )
    positions = [text.index(heading) for heading in headings]
    assert positions == sorted(positions)
    for link in (
        "docs/wiki/getting-started/quick-start.md#4-review-qc-without-changing-the-vcf",
        "examples/configs/harmonisation/sample_sheet_quantitative.csv",
        "examples/configs/harmonisation/sample_sheet_case_control.csv",
        "docs/wiki/getting-started/migrating-from-v1.md",
    ):
        assert f"]({link})" in text
    assert "not a bundled demonstration dataset" in text
    assert "QC creates an assessment, not a filtered VCF" in text
    assert "A file's existence is not enough" in text
    for flag in (
        "--apply-filter", "--apply-imputation", "--apply-manhattan", "--heritability",
    ):
        assert f"| `{flag}` |" in text


def test_readme_sample_sheet_summary_preserves_source_requirements():
    text = " ".join(ROOT_README.read_text(encoding="utf-8").split())
    assert (
        "python -m postgwas.modules.harmonisation.sample_sheet_generator"
        " \\ --input-directory raw_data \\ --output studies.csv"
    ) in text
    assert "Successful draft generation does not mean the studies are ready to run" in text
    assert "Exactly one internal `effect_allele_frequency_column` or external file/column pair" in text
    assert "Internal INFO takes priority over an external INFO file/column pair" in text
    assert "With neither source, explicitly choose `--fixed-info VALUE`" in text
    assert "`control_count_column` or `control_count` represents total N" in text
    assert "provide real control and case counts" in text


def test_readme_dependency_table_matches_method_aware_planner():
    text = ROOT_README.read_text(encoding="utf-8")
    section = text.split("### Pipeline targets and execution order", 1)[1].split(
        "### Optional pipeline stages", 1,
    )[0]
    rows = {}
    for line in section.splitlines():
        if line.startswith("| `"):
            label, plan = [cell.strip() for cell in line.strip("|").split("|")]
            rows[label] = plan.strip("`")
    public_targets = {
        name for name in REGISTRY.names()
        if REGISTRY.get(name).pipeline_enabled and REGISTRY.get(name).runner
    }
    documented_targets = {re.findall(r"`([^`]+)`", label)[0] for label in rows}
    assert documented_targets == public_targets
    guide = (REPOSITORY_ROOT / "docs/wiki/core/pipeline-workflow.md").read_text(
        encoding="utf-8",
    )
    for label, plan in rows.items():
        assert f"| {label} | `{plan}` |" in guide

    # Select the same methods named in the tables, not the packaged defaults.
    cases = [
        ("`ld_clump`, `standard` only", "ld_clump", ("--clumping-methods", "standard")),
        ("`ld_clump`, `region`", "ld_clump", ("--clumping-methods", "region")),
        ("`ld_clump`, `cojo-slct`", "ld_clump", ("--clumping-methods", "cojo-slct")),
        ("`single_cell`, `magma_celltype` or `scdrs`", "single_cell", ("--tools", "magma_celltype")),
        ("`single_cell`, `magma_celltype` or `scdrs`", "single_cell", ("--tools", "scdrs")),
        ("`single_cell`, `ldsc_celltype` only", "single_cell", ("--tools", "ldsc_celltype")),
    ]
    cases.extend(
        (
            f"`{target}`, standard clumping", target,
            ("--clumping-methods", "standard", "--finemap-method", engine),
        )
        for target in ("finemap", "caldera", "flames")
        for engine in ("susie", "finemap")
    )
    method_dependent_targets = {target for _label, target, _arguments in cases}
    cases.extend(
        (f"`{target}`", target, ())
        for target in sorted(public_targets - method_dependent_targets)
    )
    assert {label for label, _target, _arguments in cases} == set(rows)
    configuration = load_configuration()
    for label, target, arguments in cases:
        overrides = resolve_pipeline_dependency_overrides(arguments, configuration)
        if target == "single_cell":
            # The CLI resolves --tools with this shared helper after the factories.
            overrides[target] = single_cell_pipeline_dependencies(arguments[1:])
        plan = build_pipeline_plan([target], dependency_overrides=overrides)
        assert rows[label] == " → ".join(plan.steps)


def test_migration_guide_distinguishes_legacy_examples_from_current_contracts():
    guide = (
        REPOSITORY_ROOT / "docs/wiki/getting-started/migrating-from-v1.md"
    ).read_text(encoding="utf-8")
    text = " ".join(guide.split())
    for required in (
        "https://github.com/JIBINJOHNV/postgwas#readme",
        "`--config` for the harmonisation CSV | `--sample-sheet`",
        "`--defaults` for the harmonisation YAML | `--run-config`",
        "`config_version` equal to `2`",
        "`external_eaf_file` and `external_eaf_column`",
        "explicit `{chromosome}` path template",
        "The legacy `jibinjv/postgwas:1.3` image is not the current v2 build",
        "file existence alone is insufficient",
        "do not add headers manually",
        "does not establish numerical equivalence",
    ):
        assert required in text, required


def test_harmonisation_policy_guide_preserves_decisions_and_resource_contracts():
    text = HARMONISATION_POLICY_SOURCE.read_text(encoding="utf-8")
    for heading in (
        "## Required source choices",
        "## Module-level YAML choices outside `policies`",
        "## Complete harmonisation policy registry",
    ):
        assert heading in text
    for decision in (
        "reference_unmatched_retained_reverse_complement",
        "palindromic_ambiguous",
        "Neff = 4 / (1/Ncase + 1/Ncontrol)",
        "pvalue.zero_missing_se",
        "concordance_validation.failure.maximum_vcf_duplicate_records",
        "vcf_processing.liftover_tag_roles",
        "source-row",
    ):
        assert decision in text
    assert "](processing-order.md)" in text
    assert "](sample-sheet.md)" in text
    assert "](../../modules/harmonisation/configuration.md)" in text


def test_installation_docs_include_the_local_docker_build_and_runtime_contract():
    guide = REPOSITORY_ROOT / "docs/wiki/getting-started/installation.md"
    dockerfile = (REPOSITORY_ROOT / "Dockerfile").read_text(encoding="utf-8")
    for source in (ROOT_README, guide):
        text = source.read_text(encoding="utf-8")
        assert "docker build --platform linux/amd64 --load --tag postgwas:local ." in text
        assert "docker run --rm --platform linux/amd64 postgwas:local postgwas --help" in text
        assert "emulation" in text
        # Syntax validation only: never build an image or run an analysis here.
        blocks = re.findall(r"```(?:console|bash|shell)\n(.*?)```", text, re.DOTALL)
        for block in blocks:
            if "docker " in block:
                completed = subprocess.run(
                    ["bash", "-n"], input=block, capture_output=True,
                    text=True, check=False,
                )
                assert completed.returncode == 0, (source, completed.stderr)
    assert "ENTRYPOINT []" in dockerfile
    assert 'CMD ["postgwas", "--help"]' in dockerfile
    assert "RUN bash /opt/postgwas/tools/setup/verify_all_tools.sh" in dockerfile
    assert "](docs/wiki/getting-started/installation.md#docker-installation)" in (
        ROOT_README.read_text(encoding="utf-8")
    )


def test_docker_guide_explains_mounts_resources_and_validation_limits():
    guide = (REPOSITORY_ROOT / "docs/wiki/getting-started/installation.md").read_text(
        encoding="utf-8",
    )
    for required in (
        "## Docker installation",
        "dst=/input,readonly",
        "dst=/reference,readonly",
        "dst=/output",
        '--user "$(id -u):$(id -g)"',
        "--output-directory /output/qc",
        "not a published image to pull",
        "not Docker builds or full module",
        "`.gitignore` does not",
        "--mixer-backend native",
        "--mixer-backend docker",
        "MAGMA redistribution terms",
    ):
        assert required in guide, required


def test_readme_local_links_exist():
    text = ROOT_README.read_text(encoding="utf-8")
    targets = re.findall(r"\[[^\]]+\]\(([^)]+)\)", text)

    for target in targets:
        if "://" in target or target.startswith("#"):
            continue
        relative_path = target.split("#", 1)[0]
        assert (REPOSITORY_ROOT / relative_path).exists(), target


def test_documentation_workflow_validates_without_wiki_publication():
    workflow = WIKI_WORKFLOW.read_text(encoding="utf-8")

    assert workflow.count('- "README.md"') == 2
    assert workflow.count('- "tests/test_wiki_docs.py"') == 2
    assert workflow.count('- "tests/test_wiki_cli_validation.py"') == 2
    assert "name: Validate user documentation" in workflow
    assert (
        "python -m pytest -p no:cacheprovider tests/test_wiki_docs.py "
        "tests/test_wiki_cli_validation.py -q"
    ) in workflow
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
