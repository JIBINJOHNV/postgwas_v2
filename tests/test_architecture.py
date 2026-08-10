"""Regression tests for the modular orchestration boundary."""

import unittest
import subprocess
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

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
    def test_executor_uses_the_supplied_plan_without_replanning(self):
        calls = []

        def runner(args, context):
            calls.append(args._step_num)
            context["formatter"] = {"ok": True}

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused"
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))

        with patch("postgwas.pipeline.executor.REGISTRY", fake_registry), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            context = execute_pipeline(Namespace(), plan)

        self.assertEqual(calls, ["01"])
        self.assertTrue(context["formatter"]["ok"])

    def test_executor_wraps_module_failures_with_module_name(self):
        def runner(args, context):
            raise ValueError("bad input")

        fake_registry = Mock()
        fake_registry.require_pipeline_enabled.return_value = SimpleNamespace(
            runner="tests.test_architecture:unused"
        )
        plan = PipelinePlan(("formatter",), ("formatter",), ("formatter",))

        with patch("postgwas.pipeline.executor.REGISTRY", fake_registry), patch(
            "postgwas.pipeline.executor.resolve_reference", return_value=runner
        ):
            with self.assertRaisesRegex(ModuleExecutionError, "formatter: bad input"):
                execute_pipeline(Namespace(), plan)


if __name__ == "__main__":
    unittest.main()
