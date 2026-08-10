# PostGWAS modular architecture

## Engineering rules

- Keep one implementation for each shared concern and reuse it everywhere.
- Store configurable values in typed configuration; CLI options are overrides.
- Do not add backward-compatibility branches to the production package.
- Prefer small typed contracts over duplicate module-specific functions.
- Optimize I/O and parallel execution only with regression tests protecting
  scientific results and failure behavior.
- Write CLI help for users: state what to provide, when it is required, the
  default, and a short example. Keep implementation details out of help text.

## Dependency rule

Dependencies point inward:

```text
CLI -> planner/executor -> module adapters -> scientific engines
                      \-> core contracts <-/
```

`postgwas.core` must not import a scientific module. Scientific engines must not
parse CLI arguments, call `sys.exit`, or decide pipeline ordering.

Reusable implementation mechanics belong in narrowly named core modules:

- `core.dataframes` owns generic Polars casting, normalization, and counting;
- `core.values` owns scalar normalization and display formatting; and
- `core.processes` owns checked external-command execution.

A scientific module may also have a `shared/` package for plumbing that is
meaningful only inside that module. For example, `harmonisation.shared.runtime`
owns policy resolution, step contexts, logging sinks, and reject-row accounting.
Scientific formulas, evidence rules, thresholds, and interpretation stay in the
analysis file that owns them; they must not be hidden in generic utilities.

## Source layout

Scientific code lives exclusively below `postgwas.modules`. The public command
name does not have to duplicate the Python package name; the registry owns that
mapping.

```text
postgwas/
  cli/                    shared command-line components
  config/                 typed models, loaders, profiles, and YAML defaults
  core/                   contracts, errors, execution, and file I/O
  pipeline/               registry, planner, and executor
  resources/              non-configuration packaged scientific resources
  modules/
    harmonisation/
    filtering/
    formatting/exporters/
    imputation/engines/
    ld_annotation/
    ld_clumping/
    fine_mapping/engines/{finemap,susie}/
    magma/
    magmacovar/
    single_cell/
    pops/
    kpops/
    caldera/
    flames/
    ldsc/
    manhattan/
    enrichment/providers/
    mixer/
    gcta_gene/
    allele_orientation/
    qc_summary/
```

MAGMA, MAGMAcovar, single-cell integration, PoPS, and FLAMES are independent
modules rather than being forced into an inaccurate umbrella package. FINEMAP
and SuSiE are engines of the shared fine-mapping service. Alternative
imputation implementations are likewise engines of the imputation service.

GCTA fastBAT and mBAT-combo share the `gcta_gene` module because they use the
same executable, summary-statistic formatting, LD reference, validation
boundary, and result publication contract. The module supports gene-, segment-,
and custom-set fastBAT plus gene-based mBAT-combo. Selecting mBAT-combo can
expose its internally calculated fastBAT component without launching a
duplicate analysis.

Development and diagnostic scripts belong in top-level `tools/`, not in the
installable package. Module documentation belongs in `docs/modules/`.

Historical, duplicate, incomplete, and experimental files that are retained
for reference live under top-level `legacy/`. Production code, tests, package
metadata, and CLI entry points must never import from `legacy/`.

## Canonical registry

`postgwas.pipeline.registry` is the only declaration of:

- module names and descriptions;
- dependency relationships;
- standalone CLI entry points;
- pipeline runner entry points;
- parser components;
- required pipeline options;
- pipeline availability.

Entry points are stored as import strings and resolved lazily. Consequently,
`postgwas --help` and pipeline planning do not import every analysis backend.

## Planning and execution

`build_pipeline_plan()` validates every requested target and returns an immutable
`PipelinePlan`. An unknown, unavailable, cyclic, or empty request is an error.
The executor accepts only this plan and never reconstructs dependencies.

Repeated steps are explicit plan nodes. For example, a run containing imputation
and downstream analyses formats once before imputation and once after imputation.

## Module contract

New or migrated modules should have this shape:

```text
modules/<name>/
  service.py      validate() and run()
  adapters.py     external tool and file-format adapters
  cli.py          argument translation only
```

The service accepts resolved configuration, input `Artifact` objects, and a
`RunContext`. It returns one `ModuleResult`; it does not mutate an
`argparse.Namespace`. `Artifact.metadata` must carry relevant scientific
compatibility fields such as genome build, ancestry, sample identifier, schema
version, and effect-allele convention. `RunContext` is the single shared result
store for standalone and pipeline execution.

## Configuration boundary

Typed configuration is owned by `postgwas.config`. It uses one precedence rule:

```text
CLI override > run YAML > scientific profile > packaged defaults
```

Per-module values live in `config/defaults/modules/*.yaml`; their strict models
live in `config/models/modules/*.py`. Modules must not maintain another default
registry. All resolved values are written to the run manifest. Algorithmic invariants stay
in code; scientific policies, resource locations, external commands, thresholds,
retry behavior, output controls, and resource limits belong to configuration.

## Migration acceptance criteria

A migrated module is complete only when:

1. standalone and one-step pipeline modes call the same service;
2. inputs and outputs use artifact contracts;
3. no library function calls `print()` or `sys.exit()`;
4. external commands use the shared command executor;
5. every failure reaches the shared logger and run manifest;
6. expected outputs are validated before success is returned;
7. contract, failure, and scientific regression tests pass.
