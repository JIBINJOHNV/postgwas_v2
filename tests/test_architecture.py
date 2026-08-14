"""Regression tests for the modular orchestration boundary."""

import unittest
import subprocess
import sys
from argparse import Namespace
from tempfile import TemporaryDirectory
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from rich.console import Console

from postgwas.config import load_configuration
from postgwas.core.contracts import Artifact, ModuleResult, RunContext
from postgwas.core.errors import ModuleExecutionError, PipelinePlanningError
from postgwas.pipeline.executor import execute_pipeline
from postgwas.pipeline.planner import PipelinePlan, build_pipeline_plan
from postgwas.pipeline.registry import REGISTRY


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


class RegistryTests(unittest.TestCase):
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

    def test_kpops_owns_its_detailed_direct_progress(self):
        direct_progress_owners = {
            command
            for command, spec in REGISTRY.commands().items()
            if spec.owns_direct_progress
        }
        self.assertEqual(direct_progress_owners, {"kpops"})

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
