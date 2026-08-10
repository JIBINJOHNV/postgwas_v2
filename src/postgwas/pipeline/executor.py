"""Execute a validated pipeline plan.

Planning is deliberately absent from this module.  The CLI and programmatic
callers must provide the same immutable ``PipelinePlan``, preventing the two
execution modes from constructing different dependency orders.
"""

from rich.console import Console

from postgwas.core.contracts import RunContext
from postgwas.core.errors import ModuleExecutionError
from postgwas.pipeline.planner import PipelinePlan
from postgwas.pipeline.registry import REGISTRY, resolve_reference


console = Console()


def execute_pipeline(args, plan: PipelinePlan) -> RunContext:
    """Execute each planned step once and return its run context.

    A module name may occur more than once in a plan when it represents a real
    data transition, such as formatting before and after imputation.
    """

    if not isinstance(plan, PipelinePlan):
        raise TypeError("execute_pipeline requires a validated PipelinePlan")

    context = RunContext()
    console.print("\n🚀 [bold green]Starting Execution Chain[/bold green]")

    for number, module_name in enumerate(plan.steps, 1):
        spec = REGISTRY.require_pipeline_enabled(module_name)
        runner = resolve_reference(spec.runner)
        console.print(f"   {number}. Running: [cyan]{module_name}[/cyan]")

        # Expose the planned step number to module output-directory builders.
        args._step_num = f"{number:02d}"
        try:
            runner(args, context)
        except Exception as exc:
            console.print(f"\n❌ [bold red]Pipeline Failed at step: {module_name}[/bold red]")
            raise ModuleExecutionError(module_name, str(exc)) from exc

    console.print("\n✅ [bold green]All tasks completed successfully.[/bold green]\n")
    return context
