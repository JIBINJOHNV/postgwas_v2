"""The requested run budget, not host RAM, constrains locus workers."""

import ast
import inspect
from argparse import Namespace

import pytest
import yaml

from postgwas.core.execution.runtime import safe_thread_count
from postgwas.modules.fine_mapping import service


@pytest.mark.parametrize("budget,expected", [(14, 1), (16, 1), (28, 2), (42, 3)])
def test_strict_workers_obey_explicit_budget(monkeypatch, budget, expected):
    monkeypatch.setattr(
        "postgwas.core.execution.runtime.auto_detect_ram_gb",
        lambda: pytest.fail("explicit budgets must not consult host memory"),
    )
    assert safe_thread_count(
        3, 14, available_ram_gb=budget, enforce_memory_budget=True, reporter=None,
    ) == expected


@pytest.mark.parametrize("budget", [None, 0, -1, float("nan"), float("inf"), 13.99])
def test_strict_workers_reject_infeasible_or_invalid_budget(budget):
    with pytest.raises(ValueError, match="--memory-gb / execution.memory_gb"):
        safe_thread_count(
            3, 14, available_ram_gb=budget,
            enforce_memory_budget=True, reporter=None,
        )


def test_non_strict_existing_callers_retain_single_worker_fallback():
    assert safe_thread_count(3, 14, available_ram_gb=10, reporter=None) == 1


@pytest.mark.parametrize("engine", ["susie", "finemap"])
def test_fine_mapping_compute_uses_yaml_and_explicit_overrides(tmp_path, engine):
    config = tmp_path / "run.yaml"
    config.write_text(yaml.safe_dump({
        "execution": {"threads": 3, "memory_gb": 16, "random_seed": 7},
        "modules": {"fine_mapping": {"engine": engine}},
    }), encoding="utf-8")
    resolved = service._resolve_fine_mapping_arguments(Namespace(run_config=str(config)))
    assert (resolved.threads, resolved.memory_gb, resolved.seed) == (3, 16, 7)
    overridden = service._resolve_fine_mapping_arguments(Namespace(
        run_config=str(config), threads=2, memory_gb=28, seed=8,
    ))
    assert (overridden.threads, overridden.memory_gb, overridden.seed) == (2, 28, 8)


@pytest.mark.parametrize("engine", ["susie", "finemap"])
def test_pipeline_budget_failure_precedes_resource_work(monkeypatch, engine):
    monkeypatch.setattr(service, "require_pipeline_input_vcf", lambda _: {})
    monkeypatch.setattr(
        service, "run_fine_mapping_resource_preflight",
        lambda _: pytest.fail("infeasible budget must stop before resource work"),
    )
    with pytest.raises(ValueError, match="cannot fit one worker"):
        service.preflight_fine_mapping_pipeline(Namespace(
            finemap_method=engine, finemap_ld_reference="unused-prefix",
            threads=3, memory_gb=1,
        ))


@pytest.mark.parametrize("engine", ["susie", "finemap"])
def test_direct_budget_failure_precedes_resource_work(engine):
    with pytest.raises(ValueError, match="cannot fit one worker"):
        service.run_fine_mapping(Namespace(
            finemap_method=engine, finemap_ld_reference="unused-prefix",
            locus_file="unused-loci", susie_input_file="unused-summary",
            finemap_in_files="unused-summary", threads=3, memory_gb=1,
        ))


@pytest.mark.parametrize("engine,expected_calls", [("susie", 2), ("finemap", 1)])
def test_all_engine_worker_planning_calls_pass_the_resolved_budget(engine, expected_calls):
    # Guard both the minimum-worker and larger dense-LD-estimate call sites.
    if engine == "susie":
        from postgwas.modules.fine_mapping.engines.susie import adapter
        function = adapter._run_parallel_susie
    else:
        from postgwas.modules.fine_mapping.engines.finemap import adapter
        function = adapter._run_finemap_pipeline
    calls = [
        node for node in ast.walk(ast.parse(inspect.getsource(function)))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id == "safe_thread_count"
    ]
    assert len(calls) == expected_calls
    for call in calls:
        keywords = {item.arg: item.value for item in call.keywords}
        assert ast.unparse(keywords["available_ram_gb"]) == "args.memory_gb"
        assert ast.literal_eval(keywords["enforce_memory_budget"]) is True
