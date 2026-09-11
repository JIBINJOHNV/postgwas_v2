"""Deterministic admission of unequal workloads under a memory budget.

Reservations are estimates supplied by the caller, not operating-system limits.
Only admitted jobs may be submitted to an executor: its internal queue must not
bypass admission. All sizes use the same caller-defined unit.
"""

from __future__ import annotations

import math
from typing import Mapping


class ResourceQueue:
    """Alternate large and small fitting jobs and refill after completion."""

    def __init__(self, costs: Mapping[str, float], budget: float, workers: int, *, mix_sizes: bool = True):
        if not math.isfinite(budget) or budget <= 0 or workers < 1:
            raise ValueError("Scheduling requires positive memory and worker budgets.")
        self.costs = dict(costs)
        for name, cost in self.costs.items():
            if not math.isfinite(cost) or cost <= 0:
                raise ValueError("Invalid memory estimate for %s: %r" % (name, cost))
            if cost > budget:
                raise ValueError(
                    "%s requires an estimated %.3f GB, exceeding the %.3f GB "
                    "worker budget. Increase --memory-gb or review the configured "
                    "memory estimates before starting analysis." % (name, cost, budget)
                )
        self.budget = budget
        self.workers = workers
        self.pending = (
            sorted(self.costs, key=lambda k: (-self.costs[k], k))
            if mix_sizes else list(self.costs)
        )
        self.mix_sizes = mix_sizes
        self.active: dict[str, float] = {}
        self.prefer_large = True

    @property
    def reserved(self) -> float:
        return math.fsum(self.active.values())

    def admit(self) -> str | None:
        """Reserve one fitting job; None means wait for a running job."""
        if len(self.active) >= self.workers:
            return None
        order = self.pending if self.prefer_large else reversed(self.pending)
        for name in order:
            if math.fsum([*self.active.values(), self.costs[name]]) <= self.budget:
                self.pending.remove(name)
                self.active[name] = self.costs[name]
                if self.mix_sizes:
                    self.prefer_large = not self.prefer_large
                return name
        return None

    def complete(self, name: str) -> None:
        """Release a finished job, whether it succeeded or failed."""
        del self.active[name]
