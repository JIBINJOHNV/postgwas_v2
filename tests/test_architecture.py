"""Regression tests for the modular orchestration boundary."""

import os
import subprocess
import sys
import unittest
from argparse import Namespace
from tempfile import TemporaryDirectory
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pytest
from rich.console import Console
import yaml

from postgwas.config import load_configuration
from postgwas.core.contracts import Artifact, ModuleResult, RunContext
from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.errors import (
    MissingRequiredArgumentsError,
    ModuleExecutionError,
    PipelinePlanningError,
)
from postgwas.pipeline.executor import execute_pipeline
from postgwas.pipeline.cli import (
    _print_pipeline_preflight_summary,
    _require_pipeline_arguments,
    _run_pipeline_preflights,
    _validate_pipeline_entry_vcf,
)
from postgwas.pipeline.planner import PipelinePlan, build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY, resolve_reference
from preflight_support import pipeline_input_vcf_evidence


class ContractTests(unittest.TestCase):
    def test_artifact_normalises_path(self):
        artifact = Artifact("gwas_vcf", "result.vcf.gz", {"genome_build": "GRCh37"})
        self.assertEqual(artifact.path.name, "result.vcf.gz")
        self.assertEqual(artifact.metadata["genome_build"], "GRCh37")

    def test_context_supports_results_and_module_values(self):
        context = RunContext({"input": {"path": "source.txt"}})
        result = ModuleResult("example")
        context.publish(result)
        self.assertEqual(context["input"]["path"], "source.txt")
        self.assertIs(context["example"], result)

    def test_context_keeps_preflight_evidence_out_of_checkpoint_snapshot(self):
        evidence = object()
        context = RunContext(
            {"input": {"path": "source.txt"}},
            validations={"formatter": evidence},
        )

        self.assertIs(context.validation("formatter"), evidence)
        self.assertEqual(context.validation_modules(), ("formatter",))
        self.assertNotIn("formatter", context.snapshot())
        context.clear()
        self.assertIs(context.validation("formatter"), evidence)

    def test_context_rejects_empty_validation_evidence(self):
        context = RunContext()
        with self.assertRaisesRegex(ValueError, "module name"):
            context.publish_validation("", object())
        with self.assertRaisesRegex(ValueError, "must not be None"):
            context.publish_validation("magma", None)


def test_imputation_service_returns_the_selected_dataset_result(monkeypatch, tmp_path):
    from postgwas.modules.imputation import service
    from postgwas.modules.imputation.service import PredLDPipelineResources

    formatted = tmp_path / "formatted"
    formatted.mkdir()
    monkeypatch.setattr(service, "run_pred_ld_parallel", lambda **_kwargs: None)
    monkeypatch.setattr(
        service,
        "process_pred_ld_results_all_parallel",
        lambda **_kwargs: (None, None, str(tmp_path / "samples.csv")),
    )
    monkeypatch.setattr(
        service,
        "run_harmonisation", lambda _args, **_kwargs: {
            "STUDY_imputed": {
                "GRCh37": "study.grch37.vcf.gz", "status": "OK",
            }
        },
    )
    configuration = load_configuration(cli_overrides={
        "run.dataset_id": "STUDY",
        "run.output_directory": str(tmp_path / "imputation"),
        "resources.root": str(tmp_path),
        "modules.imputation.genome_build": "GRCh37",
        "modules.imputation.ld_reference_directory": str(tmp_path / "reference"),
    })
    resources = PredLDPipelineResources(
        configuration=configuration,
        reference_directory=tmp_path / "reference",
        reference_files=(),
        pred_ld_script=Path(__file__),
        harmonisation_executables={"python": sys.executable},
        resource_root=tmp_path,
        file_identities=(),
    )
    args = Namespace(
        imputation_engine="pred_ld",
        pred_ld_input_directory=str(formatted),
        imputation_ld_reference=str(tmp_path / "reference"),
        output_directory=str(tmp_path / "imputation"),
        dataset_id="STUDY",
        imputation_r2_threshold=0.8,
        imputation_minimum_maf=0.001,
        population="EUR",
        ref="TOP_LD",
        corr_method="pearson",
        threads=1,
        memory_gb=1,
    )

    result = service.run_sumstat_imputation_direct(
        args, pipeline_resources=resources,
    )

    assert result == {"GRCh37": "study.grch37.vcf.gz", "status": "OK"}


def test_imputation_runner_selects_the_validated_input_build(monkeypatch, tmp_path):
    from postgwas.pipeline import runners
    from postgwas.modules.imputation.service import PredLDPipelineResources

    resources = PredLDPipelineResources(
        configuration=load_configuration(cli_overrides={
            "resources.root": str(tmp_path),
            "modules.imputation.genome_build": "GRCh38",
            "modules.imputation.ld_reference_directory": str(tmp_path / "reference"),
        }),
        reference_directory=tmp_path / "reference",
        reference_files=(),
        pred_ld_script=Path(__file__),
        harmonisation_executables={"python": sys.executable},
        resource_root=tmp_path,
        file_identities=(),
    )

    context = RunContext(
        {"formatter": {"pred_ld": {"pred_ld_folder": "formatted"}}},
        validations={
            "current_vcf": {
                "indexed": SimpleNamespace(genome_build="GRCh38"),
            },
            "imputation": PipelinePreflightEvidence(
                "imputation", {}, resources,
            ),
        },
    )
    args = Namespace(
        output_directory=str(tmp_path),
        dataset_id="STUDY",
        _step_num=2,
    )
    outputs = {
        "GRCh37": "study.grch37.vcf.gz",
        "GRCh38": "study.grch38.vcf.gz",
    }
    validated = []
    monkeypatch.setattr(
        "postgwas.modules.imputation.service.run_sumstat_imputation_direct",
        lambda _args, **_kwargs: outputs,
    )
    monkeypatch.setattr(
        runners,
        "_validate_current_pipeline_vcf",
        lambda current_args, _ctx: validated.append(current_args.vcf),
    )

    result = runners.run_imputation_runner(args, context)

    assert result == outputs
    assert args.vcf == "study.grch38.vcf.gz"
    assert validated == ["study.grch38.vcf.gz"]


def test_current_vcf_validation_reuses_the_entry_bcftools(monkeypatch):
    from postgwas.pipeline import runners

    entry = pipeline_input_vcf_evidence()["input_vcf"]
    context = RunContext(validations={
        "input_vcf": entry,
        "current_vcf": entry,
    })
    observed = {}

    def validate(vcf, dataset_id, bcftools, _module, *, cached):
        observed.update({
            "vcf": vcf,
            "dataset_id": dataset_id,
            "bcftools": bcftools,
            "cached": cached,
        })
        return entry

    monkeypatch.setattr(
        "postgwas.core.vcf.validate_harmonised_indexed_vcf",
        validate,
    )
    monkeypatch.setattr(
        runners,
        "resolve_executable",
        lambda *_args, **_kwargs: pytest.fail(
            "bcftools was resolved again after entry preflight"
        ),
    )

    result = runners._validate_current_pipeline_vcf(
        Namespace(vcf="study.vcf.gz", dataset_id="STUDY"), context,
    )

    assert result is entry
    assert observed["bcftools"] == "/validated/bcftools"
    assert observed["cached"] is entry["indexed"]
    assert context.validation("current_vcf") is entry


def test_entry_vcf_validation_retains_bcftools_identity(monkeypatch, tmp_path):
    executable = tmp_path / "bcftools"
    executable.write_text("validated executable\n", encoding="utf-8")
    expected = pipeline_input_vcf_evidence()["input_vcf"]
    monkeypatch.setattr(
        "postgwas.pipeline.cli.resolve_executable",
        lambda *_args, **_kwargs: str(executable),
    )
    monkeypatch.setattr(
        "postgwas.core.vcf.validate_harmonised_indexed_vcf",
        lambda *_args, **_kwargs: dict(expected),
    )

    evidence = _validate_pipeline_entry_vcf(
        Namespace(vcf="study.vcf.gz", dataset_id="STUDY"),
        load_configuration(),
    )

    assert evidence["bcftools_identity"][0].path == executable.resolve()


def test_pipeline_preflight_summary_reports_ordered_readiness(monkeypatch):
    initial = pipeline_input_vcf_evidence()
    input_vcf = initial["input_vcf"]
    evidence = {
        **initial,
        "formatter": PipelinePreflightEvidence("formatter", input_vcf),
        "magma": PipelinePreflightEvidence(
            "magma",
            input_vcf,
            resources={"validated": True},
            deferred_checks=("Validate generated MAGMA input.",),
        ),
    }
    observed = []
    monkeypatch.setattr(
        "postgwas.pipeline.cli.print_screen_block", observed.append,
    )

    _print_pipeline_preflight_summary(
        input_vcf, evidence, ("formatter", "magma"),
    )

    screen = observed[0]
    assert "Pipeline readiness" in screen
    assert "1 variants; indexed single-sample VCF validated" in screen
    assert "Stage preflights" in screen and "2/2 passed" in screen
    assert "1 resource-dependent stages passed" in screen
    assert "1 deferred to their consuming stages" in screen


def test_kpops_runner_forwards_resource_preflight(monkeypatch, tmp_path):
    from postgwas.pipeline import runners

    resources = object()
    context = RunContext(
        {
            "magma": {
                "primary_mapping": "positional",
                "mapping_analyses": {
                    "positional": {
                        "result_statistic_type": "calibrated_gene_p_value",
                    },
                },
                "magma_genes_prefix": "validated_magma",
            },
        },
        validations={
            "kpops": PipelinePreflightEvidence("kpops", {}, resources),
        },
    )
    observed = {}

    def run(args, ctx, *, pipeline_resources):
        observed.update({
            "prefix": args.magma_association_prefix,
            "context": ctx,
            "resources": pipeline_resources,
        })
        return {"status": "success"}

    monkeypatch.setattr(
        "postgwas.modules.kpops.service.run_kpops_direct", run,
    )
    args = Namespace(output_directory=str(tmp_path), _step_num=1)

    result = runners.run_kpops_runner(args, context)

    assert result == {"status": "success"}
    assert observed == {
        "prefix": "validated_magma",
        "context": context,
        "resources": resources,
    }


def test_caldera_runner_forwards_resource_preflight(monkeypatch, tmp_path):
    from postgwas.pipeline import runners

    resources = object()
    context = RunContext(
        {
            "pops_output": "validated_pops.preds",
            "finemap": {"flames_input": "validated_credible_sets"},
        },
        validations={
            "caldera": PipelinePreflightEvidence("caldera", {}, resources),
        },
    )
    observed = {}

    def run(args, ctx, *, pipeline_resources):
        observed.update({
            "pops": args.pops_file,
            "credible_sets": args.finemap_credible_sets_directory,
            "context": ctx,
            "resources": pipeline_resources,
        })
        return {"status": "success"}

    monkeypatch.setattr(
        "postgwas.modules.caldera.service.run_caldera_direct", run,
    )
    args = Namespace(output_directory=str(tmp_path), _step_num=1)

    result = runners.run_caldera_runner(args, context)

    assert result == {"status": "success"}
    assert observed == {
        "pops": "validated_pops.preds",
        "credible_sets": "validated_credible_sets",
        "context": context,
        "resources": resources,
    }


class RegistryTests(unittest.TestCase):
    def test_every_executable_pipeline_step_has_shared_preflight(self):
        missing = [
            name
            for name in REGISTRY.names(include_internal=True)
            if REGISTRY.get(name).pipeline_enabled
            and REGISTRY.get(name).runner is not None
            and REGISTRY.get(name).preflight is None
        ]
        self.assertEqual(missing, [])

    def test_every_pipeline_preflight_requires_common_vcf_evidence_first(self):
        for name in REGISTRY.names(include_internal=True):
            specification = REGISTRY.get(name)
            if (
                not specification.pipeline_enabled
                or specification.runner is None
                or specification.preflight is None
            ):
                continue
            validator = resolve_reference(specification.preflight)
            with self.subTest(module=name), self.assertRaisesRegex(
                RuntimeError,
                "shared input-VCF validation",
            ):
                validator(Namespace(), preflight_evidence={})

    def test_vcf_only_pipeline_modules_have_shared_preflights(self):
        expected = {
            "sumstat_filter": "preflight_sumstat_filter",
            "post_imputation_filter": "preflight_post_imputation_filter",
            "formatter": "preflight_formatter",
            "manhattan": "preflight_manhattan",
            "qc_summary": "preflight_qc_summary",
        }
        for module, function in expected.items():
            self.assertEqual(
                REGISTRY.get(module).preflight,
                "postgwas.pipeline.preflight:%s" % function,
            )

    def test_pipeline_preflights_return_evidence_in_module_order(self):
        initial = pipeline_input_vcf_evidence()
        validators = {
            "first:validate": lambda _args, **_kwargs: PipelinePreflightEvidence(
                "first", initial["input_vcf"], {"validated": "first"},
            ),
            "second:validate": lambda _args, **_kwargs: PipelinePreflightEvidence(
                "second", initial["input_vcf"], {"validated": "second"},
            ),
        }

        def specification(name):
            return SimpleNamespace(preflight="%s:validate" % name)

        with patch(
            "postgwas.pipeline.cli.REGISTRY.get", side_effect=specification,
        ), patch(
            "postgwas.pipeline.cli.resolve_reference",
            side_effect=lambda reference: validators[reference],
        ):
            evidence = _run_pipeline_preflights(
                Namespace(), ("first", "second"), initial_evidence=initial,
            )

        self.assertEqual(tuple(evidence), ("input_vcf", "first", "second"))
        self.assertEqual(
            evidence["first"].resources, {"validated": "first"},
        )

    def test_pipeline_rejects_a_preflight_that_discards_its_evidence(self):
        with patch(
            "postgwas.pipeline.cli.REGISTRY.get",
            return_value=SimpleNamespace(preflight="example:validate"),
        ), patch(
            "postgwas.pipeline.cli.resolve_reference",
            return_value=lambda _args, **_kwargs: None,
        ):
            with self.assertRaisesRegex(
                RuntimeError, "PipelinePreflightEvidence contract",
            ):
                _run_pipeline_preflights(Namespace(), ("example",))

    def test_pipeline_required_arguments_use_shared_resolved_contract(self):
        configuration = load_configuration()
        arguments = Namespace(
            vcf=None,
            ld_folder=None,
            finemap_ld_reference=None,
        )

        with self.assertRaisesRegex(
            MissingRequiredArgumentsError,
            r"Required argument not provided: --vcf\. Provide --vcf VALUE\.",
        ) as captured:
            _require_pipeline_arguments(
                arguments,
                ("ld_clump", "formatter", "finemap"),
                configuration,
                explicitly_requested_modules=("finemap",),
            )

        message = str(captured.exception)
        self.assertEqual(message.count("--vcf. Provide --vcf VALUE."), 1)
        self.assertIn("modules.ld_clumping.reference.directory", message)
        self.assertIn("modules.fine_mapping.input.ld_reference_prefix", message)

    def test_top_level_dispatcher_is_lazy(self):
        code = (
            "import sys; import postgwas.__main__; "
            "assert 'postgwas.modules.flames.service' not in sys.modules; "
            "assert 'postgwas.modules.imputation.engines.pred_ld.pred_ld_runner' "
            "not in sys.modules"
        )
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_pipeline_runners_do_not_import_unselected_backends(self):
        code = (
            "import sys; import postgwas.pipeline.runners; "
            "assert 'postgwas.modules.magma.service' not in sys.modules; "
            "assert 'statsmodels' not in sys.modules"
        )
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_mixer_is_available_standalone_and_in_pipeline(self):
        commands = REGISTRY.commands()
        self.assertEqual(
            commands["mixer"].cli_entrypoint,
            "postgwas.modules.mixer.cli:main",
        )
        self.assertTrue(commands["mixer"].pipeline_enabled)

    def test_every_enabled_module_has_a_runner(self):
        for name in REGISTRY.names(include_internal=True):
            spec = REGISTRY.get(name)
            if spec.pipeline_enabled:
                self.assertTrue(spec.runner, name)

    def test_scientific_entrypoints_use_the_modules_namespace(self):
        for command, spec in REGISTRY.commands().items():
            if command in {"config", "pipeline", "resources"}:
                continue
            self.assertTrue(
                spec.cli_entrypoint.startswith("postgwas.modules."),
                (command, spec.cli_entrypoint),
            )

    def test_only_pipeline_executor_owns_its_top_level_progress(self):
        owners = {
            command
            for command, spec in REGISTRY.commands().items()
            if spec.owns_top_level_progress
        }
        self.assertEqual(owners, {"pipeline"})

    def test_modules_with_detailed_direct_stages_own_their_progress(self):
        direct_progress_owners = {
            command
            for command, spec in REGISTRY.commands().items()
            if spec.owns_direct_progress
        }
        self.assertEqual(direct_progress_owners, {"caldera", "kpops", "pops"})

    def test_every_scientific_direct_command_has_a_checkpoint_strategy(self):
        strategies = {
            command: spec.direct_checkpoint
            for command, spec in REGISTRY.commands().items()
        }
        self.assertEqual(
            {
                command for command, strategy in strategies.items()
                if strategy == "not_applicable"
            },
            {"config", "pipeline", "resources"},
        )
        self.assertTrue(all(
            strategy in {"orchestrated", "native", "not_applicable"}
            for strategy in strategies.values()
        ))
        self.assertEqual(strategies["formatter"], "orchestrated")
        self.assertEqual(strategies["magma"], "native")

    def test_unknown_module_fails_during_planning(self):
        with self.assertRaisesRegex(PipelinePlanningError, "Unknown module 'unknown'"):
            build_pipeline_plan(["unknown"])

    def test_mixer_plan_includes_formatter(self):
        self.assertEqual(build_pipeline_plan(["mixer"]).steps, ("formatter", "mixer"))

    def test_gcta_gene_plan_includes_formatter(self):
        self.assertEqual(
            build_pipeline_plan(["gcta_gene"]).steps,
            ("formatter", "gcta_gene"),
        )

    def test_gcta_cojo_plan_includes_formatter(self):
        self.assertEqual(
            build_pipeline_plan(["gcta_cojo"]).steps,
            ("formatter", "gcta_cojo"),
        )

    def test_single_cell_plan_reuses_formatter_and_magma(self):
        self.assertEqual(
            build_pipeline_plan(["single_cell"]).steps,
            ("formatter", "magma", "single_cell"),
        )


class LayoutTests(unittest.TestCase):
    def test_canonical_module_packages_exist(self):
        package_root = Path(__file__).parents[1] / "src" / "postgwas" / "modules"
        expected = {
            "allele_orientation",
            "caldera",
            "enrichment",
            "filtering",
            "fine_mapping",
            "flames",
            "formatting",
            "gcta_cojo",
            "gcta_gene",
            "harmonisation",
            "imputation",
            "kpops",
            "ld_annotation",
            "ld_clumping",
            "ldsc",
            "magma",
            "magmacovar",
            "manhattan",
            "mixer",
            "pops",
            "qc_summary",
            "single_cell",
        }
        discovered = {
            path.name
            for path in package_root.iterdir()
            if path.is_dir() and (path / "__init__.py").exists()
        }
        self.assertEqual(discovered, expected)

    def test_multi_tool_modules_use_engine_packages(self):
        package_root = Path(__file__).parents[1] / "src" / "postgwas" / "modules"
        self.assertTrue(
            (package_root / "fine_mapping" / "engines" / "finemap").is_dir()
        )
        self.assertTrue(
            (package_root / "fine_mapping" / "engines" / "susie").is_dir()
        )
        self.assertTrue((package_root / "imputation" / "engines").is_dir())

    def test_production_package_contains_no_legacy_source_names(self):
        source_root = Path(__file__).parents[1] / "src" / "postgwas"
        forbidden_parts = {"clis", "scripts", "utils"}
        forbidden_fragments = (
            "_back.py",
            "_old.py",
            "_original.py",
            "cli_withpipeline.py",
            "cli_withpipelines.py",
            "need_to_delete",
        )

        offenders = []
        for path in source_root.rglob("*.py"):
            relative = path.relative_to(source_root)
            if forbidden_parts.intersection(relative.parts):
                offenders.append(str(relative))
                continue
            if any(fragment in path.name.lower() for fragment in forbidden_fragments):
                offenders.append(str(relative))

        self.assertEqual(offenders, [])


class PlannerTests(unittest.TestCase):
    def test_empty_target_is_invalid(self):
        with self.assertRaisesRegex(PipelinePlanningError, "at least one"):
            build_pipeline_plan([])

    def test_explicit_formatter_produces_one_step(self):
        plan = build_pipeline_plan(["formatter"])
        self.assertEqual(plan.steps, ("formatter",))

    def test_finemap_dependencies_are_in_execution_order(self):
        plan = build_pipeline_plan(["finemap"])
        self.assertEqual(
            plan.steps,
            ("annot_ldblock", "ld_clump", "formatter", "finemap"),
        )

    def test_standard_finemap_omits_ld_annotation_and_clumps_first(self):
        from postgwas.pipeline.planner import (
            resolve_pipeline_dependency_overrides,
        )

        overrides = resolve_pipeline_dependency_overrides(
            ("--clumping-methods", "standard"),
            load_configuration(),
        )
        self.assertEqual(overrides["ld_clump"], ())
        self.assertEqual(
            build_pipeline_plan(
                ["finemap"], dependency_overrides=overrides,
            ).steps,
            ("ld_clump", "formatter", "finemap"),
        )

    def test_region_and_cojo_add_only_their_real_clumping_dependencies(self):
        from postgwas.pipeline.planner import (
            resolve_pipeline_dependency_overrides,
        )

        configuration = load_configuration()
        region = resolve_pipeline_dependency_overrides(
            ("--clumping-methods", "region", "standard"),
            configuration,
        )
        cojo = resolve_pipeline_dependency_overrides(
            ("--clumping-methods", "standard", "cojo-slct"),
            configuration,
        )
        combined = resolve_pipeline_dependency_overrides(
            (
                "--clumping-methods", "region", "standard", "cojo-slct",
            ),
            configuration,
        )

        self.assertEqual(region["ld_clump"], ("annot_ldblock",))
        self.assertEqual(cojo["ld_clump"], ("formatter",))
        self.assertEqual(
            build_pipeline_plan(
                ["finemap"], dependency_overrides=region,
            ).steps,
            ("annot_ldblock", "ld_clump", "formatter", "finemap"),
        )
        self.assertEqual(
            build_pipeline_plan(
                ["finemap"], dependency_overrides=cojo,
            ).steps,
            ("formatter", "ld_clump", "finemap"),
        )
        self.assertEqual(
            build_pipeline_plan(
                ["finemap"], dependency_overrides=combined,
            ).steps,
            ("annot_ldblock", "formatter", "ld_clump", "finemap"),
        )

    def test_filter_and_imputation_have_pre_and_post_filter_stages(self):
        plan = build_pipeline_plan([], apply_filter=True, apply_imputation=True)
        self.assertEqual(
            plan.steps,
            ("sumstat_filter", "formatter", "imputation", "post_imputation_filter"),
        )

    def test_imputation_and_downstream_analysis_format_twice(self):
        plan = build_pipeline_plan(["magma"], apply_imputation=True)
        self.assertEqual(plan.steps.count("formatter"), 2)
        self.assertLess(plan.steps.index("imputation"), plan.steps.index("magma"))

    def test_imputation_formatter_stages_have_distinct_consumers(self):
        from postgwas.modules.formatting.contracts import required_formats
        from postgwas.pipeline.runners import (
            _formatter_consumers_for_current_vcf,
        )

        args = Namespace(
            modules=["formatter", "imputation", "magma", "mixer"],
            apply_imputation=True,
        )

        before, is_pre_imputation = _formatter_consumers_for_current_vcf(
            args, RunContext(),
        )
        after, is_post_imputation = _formatter_consumers_for_current_vcf(
            args, RunContext({"imputation": {"GRCh37": "imputed.vcf.gz"}}),
        )

        self.assertEqual(before, ["imputation"])
        self.assertTrue(is_pre_imputation)
        self.assertEqual(after, ["formatter", "magma", "mixer"])
        self.assertFalse(is_post_imputation)

        formatting = load_configuration().modules.formatting
        self.assertEqual(required_formats(formatting, before), ["pred_ld"])
        self.assertEqual(
            required_formats(formatting, after), ["magma", "mixer"],
        )

    def test_plan_is_immutable(self):
        plan = build_pipeline_plan(["formatter"])
        with self.assertRaises(AttributeError):
            plan.steps = ()


class ExecutorTests(unittest.TestCase):
    @staticmethod
    def _configuration(output_directory):
        return load_configuration(cli_overrides={
            "run.output_directory": str(output_directory),
        })

    def test_executor_uses_the_supplied_plan_without_replanning(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("ok\n", encoding="utf-8")
            context["formatter"] = {"ok": True, "output": str(output)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        stream = StringIO()
        progress_console = Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        )

        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ), patch("postgwas.pipeline.executor.console", progress_console):
            args = Namespace(
                output_directory=directory, resume=True, overwrite=False,
            )
            context = execute_pipeline(
                args, plan, self._configuration(directory),
            )

        self.assertEqual(calls, ["01"])
        self.assertTrue(context["formatter"]["ok"])
        self.assertIn("Completed 1/1 · Convert GWAS-VCF inputs", stream.getvalue())

    def test_output_directory_name_is_not_discovered_as_a_stage_input(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            stage = Path(args.output_directory) / "01_formatter"
            log = stage / "logs" / "STUDY_formatter.log"
            with log.open("a", encoding="utf-8") as handle:
                handle.write("current run\n")
            output = stage / "result.tsv"
            output.write_text("validated\n", encoding="utf-8")
            context["formatter"] = {
                "output": str(output),
                "log_file": str(log),
            }
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(
            ("magmacovar",),
            ("formatter", "magma", "magmacovar"),
            ("formatter",),
        )
        stream = StringIO()
        progress_console = Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        )

        with TemporaryDirectory() as directory:
            working_directory = Path(directory)
            output_directory = working_directory / "magmacovar"
            log = output_directory / "01_formatter/logs/STUDY_formatter.log"
            log.parent.mkdir(parents=True)
            log.write_text("prior run\n", encoding="utf-8")
            original_directory = Path.cwd()
            try:
                os.chdir(working_directory)
                with patch(
                    "postgwas.pipeline.executor.REGISTRY", fake_registry,
                ), patch(
                    "postgwas.pipeline.executor.resolve_reference",
                    return_value=runner,
                ), patch(
                    "postgwas.pipeline.executor.console", progress_console,
                ):
                    context = execute_pipeline(
                        Namespace(
                            output_directory=str(output_directory),
                            modules=["magmacovar"],
                            resume=True,
                            overwrite=False,
                        ),
                        plan,
                        self._configuration(output_directory),
                    )
                checkpoint_document = yaml.safe_load(
                    (
                        output_directory
                        / "run_metadata/checkpoints/01_formatter.yaml"
                    ).read_text(encoding="utf-8")
                )
            finally:
                os.chdir(original_directory)

        self.assertEqual(calls, ["01"])
        self.assertEqual(
            context["formatter"]["output"],
            str(output_directory / "01_formatter/result.tsv"),
        )
        self.assertNotIn(str(log.resolve()), checkpoint_document["inputs"])
        self.assertIn(
            str(log.resolve()),
            {
                record["path"]
                for record in checkpoint_document["outputs"].values()
            },
        )
        self.assertIn(
            "Completed 1/1 · Convert GWAS-VCF inputs",
            stream.getvalue(),
        )

    def test_live_screen_log_below_module_directory_is_not_a_stage_input(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            with screen_log.open("a", encoding="utf-8") as handle:
                handle.write("formatter progress\n")
            stage = output_directory / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            result = stage / "result.tsv"
            result.write_text("validated\n", encoding="utf-8")
            context["formatter"] = {"output": str(result)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
            pipeline_progress_factory=None,
            pipeline_title_factory=None,
        )
        plan = PipelinePlan(
            ("gcta_gene",),
            ("formatter", "gcta_gene"),
            ("formatter",),
        )
        stream = StringIO()
        progress_console = Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        )

        with TemporaryDirectory() as directory:
            study_root = Path(directory) / "SCZ_2026"
            output_directory = study_root / "gcta_gene/mbat_combo"
            screen_log = output_directory / "run_metadata/screen.log"
            screen_log.parent.mkdir(parents=True)
            screen_log.write_text("pipeline start\n", encoding="utf-8")
            original_directory = Path.cwd()
            try:
                os.chdir(study_root)
                with patch(
                    "postgwas.pipeline.executor.REGISTRY", fake_registry,
                ), patch(
                    "postgwas.pipeline.executor.resolve_reference",
                    return_value=runner,
                ), patch(
                    "postgwas.pipeline.executor.console", progress_console,
                ):
                    execute_pipeline(
                        Namespace(
                            output_directory=str(output_directory),
                            modules=["formatter", "gcta_gene"],
                            method="mbat_combo",
                            resume=True,
                            overwrite=False,
                        ),
                        plan,
                        self._configuration(output_directory),
                    )
                checkpoint_document = yaml.safe_load(
                    (
                        output_directory
                        / "run_metadata/checkpoints/01_formatter.yaml"
                    ).read_text(encoding="utf-8")
                )
            finally:
                os.chdir(original_directory)

        self.assertEqual(calls, ["01"])
        self.assertNotIn(
            str(screen_log.resolve()), checkpoint_document["inputs"],
        )
        self.assertIn("Completed 1/1", stream.getvalue())

    def test_executor_wraps_module_failures_with_module_name(self):
        def runner(args, context):
            raise ValueError("bad input")

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        stream = StringIO()
        progress_console = Console(
            file=stream,
            force_terminal=True,
            color_system=None,
            width=120,
        )

        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ), patch("postgwas.pipeline.executor.console", progress_console):
            args = Namespace(
                output_directory=directory, resume=True, overwrite=False,
            )
            with self.assertRaisesRegex(ModuleExecutionError, "formatter: bad input"):
                execute_pipeline(
                    args, plan, self._configuration(directory),
                )
        self.assertIn("Failed 1/1 · Convert GWAS-VCF inputs", stream.getvalue())
        self.assertNotIn("All 1 stages completed", stream.getvalue())

    def test_executor_uses_one_detailed_progress_plan_across_pipeline_modules(self):
        stages = tuple("Scientific stage %d" % number for number in range(1, 13))
        calls = []

        def factory(_args):
            return {
                "label": "Detailed scientific progress",
                "stages": stages,
                "modules": {"formatter": (1, 6), "magma": (6, 12)},
                "deferred_completion_modules": {"formatter": 6},
            }

        def formatter(args, context):
            calls.append("formatter")
            for number in range(1, 6):
                args._pipeline_stage_progress.start(number)
                args._pipeline_stage_progress.complete(number)
            args._pipeline_stage_progress.start(6)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "formatter.tsv"
            output.write_text("ok\n", encoding="utf-8")
            context["formatter"] = {"output": str(output)}
            return context["formatter"]

        def magma(args, context):
            calls.append("magma")
            args._pipeline_stage_progress.complete(6)
            for number in range(7, 13):
                args._pipeline_stage_progress.start(number)
                args._pipeline_stage_progress.complete(number)
            stage = Path(args.output_directory) / "02_magma"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "magma.tsv"
            output.write_text("ok\n", encoding="utf-8")
            context["magma"] = {"output": str(output)}
            return context["magma"]

        specs = {
            "formatter": SimpleNamespace(
                runner="formatter", description="formatter",
                pipeline_output_name=None, pipeline_title_factory=None,
            ),
            "magma": SimpleNamespace(
                runner="magma", description="magma",
                pipeline_output_name=None, pipeline_title_factory=None,
                pipeline_progress_factory="factory",
            ),
        }
        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.side_effect = specs.__getitem__
        references = {"factory": factory, "formatter": formatter, "magma": magma}
        stream = StringIO()
        progress_console = Console(
            file=stream, force_terminal=True, color_system=None, width=120,
        )
        plan = PipelinePlan(
            ("magma",), ("formatter", "magma"), ("formatter", "magma"),
        )
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry,
        ), patch(
            "postgwas.pipeline.executor.resolve_reference",
            side_effect=references.__getitem__,
        ), patch("postgwas.pipeline.executor.console", progress_console):
            for _run in range(2):
                execute_pipeline(
                    Namespace(
                        output_directory=directory, resume=True, overwrite=False,
                    ),
                    plan,
                    self._configuration(directory),
                )

        output = stream.getvalue()
        self.assertIn("Completed 12/12 · Scientific stage 12", output)
        self.assertIn("All 12 stages completed", output)
        self.assertNotIn("Current 1/2", output)
        self.assertEqual(calls, ["formatter", "magma"])

    def test_detailed_pipeline_progress_stops_at_the_failed_scientific_stage(self):
        def factory(_args):
            return {
                "label": "Detailed scientific progress",
                "stages": ("Validate input", "Validate reference", "Analyse"),
                "modules": {"formatter": (1, 2), "magma": (3, 3)},
            }

        def formatter(args, _context):
            args._pipeline_stage_progress.start(1)
            args._pipeline_stage_progress.complete(1)
            args._pipeline_stage_progress.start(2)
            raise ValueError("invalid reference")

        specs = {
            "formatter": SimpleNamespace(
                runner="formatter", description="formatter",
                pipeline_output_name=None, pipeline_title_factory=None,
            ),
            "magma": SimpleNamespace(
                runner="magma", description="magma",
                pipeline_output_name=None, pipeline_title_factory=None,
                pipeline_progress_factory="factory",
            ),
        }
        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.side_effect = specs.__getitem__
        references = {"factory": factory, "formatter": formatter}
        stream = StringIO()
        progress_console = Console(
            file=stream, force_terminal=True, color_system=None, width=120,
        )
        plan = PipelinePlan(
            ("magma",), ("formatter", "magma"), ("formatter", "magma"),
        )
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry,
        ), patch(
            "postgwas.pipeline.executor.resolve_reference",
            side_effect=references.__getitem__,
        ), patch("postgwas.pipeline.executor.console", progress_console):
            with self.assertRaisesRegex(
                ModuleExecutionError, "formatter: invalid reference",
            ):
                execute_pipeline(
                    Namespace(
                        output_directory=directory, resume=True, overwrite=False,
                    ),
                    plan,
                    self._configuration(directory),
                )

        output = stream.getvalue()
        self.assertIn("Failed 2/3 · Validate reference", output)
        self.assertNotIn("All 3 stages completed", output)
        self.assertNotIn("100%", output)

    def test_executor_resumes_completed_stage_without_rerunning(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("validated\n", encoding="utf-8")
            context["formatter"] = {"output": str(output), "rows": 1}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            configuration = self._configuration(directory)
            first = execute_pipeline(
                Namespace(
                    output_directory=directory, resume=True, overwrite=False,
                ),
                plan,
                configuration,
            )
            second = execute_pipeline(
                Namespace(
                    output_directory=directory, resume=True, overwrite=False,
                ),
                plan,
                configuration,
            )

        self.assertEqual(calls, ["01"])
        self.assertEqual(first.snapshot(), second.snapshot())

    def test_executor_changed_parameter_warns_and_restarts_stage(self):
        calls = []

        def runner(args, context):
            calls.append(args.threshold)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("threshold=%s\n" % args.threshold, encoding="utf-8")
            context["formatter"] = {"output": str(output)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            configuration = self._configuration(directory)
            for threshold in (0.05, 0.01):
                execute_pipeline(
                    Namespace(
                        output_directory=directory,
                        resume=True,
                        overwrite=False,
                        threshold=threshold,
                    ),
                    plan,
                    configuration,
                )
            audit = Path(directory) / (
                "run_metadata/checkpoints/checkpoint_events.log"
            )
            audit_text = audit.read_text(encoding="utf-8")
            result_text = (
                Path(directory) / "01_formatter" / "result.tsv"
            ).read_text(encoding="utf-8")

        self.assertEqual(calls, [0.05, 0.01])
        self.assertIn("changed_parameters", audit_text)
        self.assertEqual(result_text, "threshold=0.01\n")

    def test_executor_restarts_verified_partial_stage_after_warning(self):
        attempts = []

        def runner(args, context):
            attempts.append(len(attempts) + 1)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("attempt=%d\n" % attempts[-1], encoding="utf-8")
            if len(attempts) == 1:
                raise ValueError("simulated interruption")
            context["formatter"] = {"output": str(output)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            configuration = self._configuration(directory)
            arguments = Namespace(
                output_directory=directory, resume=True, overwrite=False,
            )
            with self.assertRaisesRegex(
                ModuleExecutionError, "simulated interruption",
            ):
                execute_pipeline(arguments, plan, configuration)
            execute_pipeline(
                Namespace(
                    output_directory=directory, resume=True, overwrite=False,
                ),
                plan,
                configuration,
            )
            audit_text = (
                Path(directory)
                / "run_metadata/checkpoints/checkpoint_events.log"
            ).read_text(encoding="utf-8")
            result_text = (
                Path(directory) / "01_formatter" / "result.tsv"
            ).read_text(encoding="utf-8")

        self.assertEqual(attempts, [1, 2])
        self.assertIn("partial_results", audit_text)
        self.assertEqual(result_text, "attempt=2\n")

    def test_overwrite_control_does_not_invalidate_the_following_resume(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            stage = Path(args.output_directory) / (
                "%s_formatter" % args._step_num
            )
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("validated\n", encoding="utf-8")
            context["formatter"] = {"output": str(output)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(
            ("formatter",),
            ("formatter",),
            ("formatter", "formatter"),
        )
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            overwrite_configuration = load_configuration(cli_overrides={
                "run.output_directory": directory,
                "run.overwrite": True,
            })
            execute_pipeline(
                Namespace(
                    output_directory=directory,
                    resume=True,
                    overwrite=True,
                ),
                plan,
                overwrite_configuration,
            )
            execute_pipeline(
                Namespace(
                    output_directory=directory,
                    resume=True,
                    overwrite=False,
                ),
                plan,
                self._configuration(directory),
            )

        self.assertEqual(calls, ["01", "02"])

    def test_executor_refuses_modified_checkpoint_output(self):
        calls = []

        def runner(args, context):
            calls.append(1)
            stage = Path(args.output_directory) / "01_formatter"
            stage.mkdir(parents=True, exist_ok=True)
            output = stage / "result.tsv"
            output.write_text("validated\n", encoding="utf-8")
            context["formatter"] = {"output": str(output)}
            return context["formatter"]

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused",
            description="Convert GWAS-VCF inputs",
            pipeline_output_name=None,
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))
        with TemporaryDirectory() as directory, patch(
            "postgwas.pipeline.executor.REGISTRY", fake_registry
        ), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            configuration = self._configuration(directory)
            execute_pipeline(
                Namespace(
                    output_directory=directory, resume=True, overwrite=False,
                ),
                plan,
                configuration,
            )
            output = Path(directory) / "01_formatter" / "result.tsv"
            output.write_text("user-modified\n", encoding="utf-8")
            with self.assertRaisesRegex(
                ModuleExecutionError, "modified after validation",
            ):
                execute_pipeline(
                    Namespace(
                        output_directory=directory,
                        resume=True,
                        overwrite=False,
                    ),
                    plan,
                    configuration,
                )
            preserved = output.read_text(encoding="utf-8")

        self.assertEqual(calls, [1])
        self.assertEqual(preserved, "user-modified\n")


if __name__ == "__main__":
    unittest.main()
