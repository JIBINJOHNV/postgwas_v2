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

Historical, duplicate, incomplete, and experimental source files are not kept
in the working tree. Git history provides recovery for superseded artifacts.
Generated analysis results, scientific reference datasets, and third-party
binaries must not be committed to the repository.

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
The first pass creates only the imputation input from the original VCF; the
second creates only the selected downstream inputs from the re-harmonised
imputed VCF. A formatter pass extracts its current VCF once and fans that table
out to every consumer assigned to that data state.

## Validation and preparation boundary

Before any pipeline stage creates scientific output, the orchestrator validates
the indexed PostGWAS harmonised GWAS-VCF once and records reusable header,
sample, build, record-count, index, and bcftools evidence. Every selected module
then runs its registered preflight against that common VCF evidence and validates
the external references, resource files, genome-build/population declarations,
and executables that are knowable before generated inputs exist. Preflights
return the shared `PipelinePreflightEvidence` contract; module-native validated
resource evidence is retained in `RunContext` for formatter and analysis
consumers.

`core.input_validation.InputValidationSession` spans startup and execution.
Registered preflights and stage runners use module scopes; shared validators
record file path, role, exact checks, status, compact metrics, explanatory
message, and consumers. Basic file checks record availability only. Existing
module evidence is adapted explicitly by `pipeline.resource_validation`; no
generic configuration walk guesses which objects are input files or which
scientific checks ran. Independent module-preflight failures are collected
together. Failure of the shared entry-VCF check blocks the dependent preflights.

Reusable file validators include PLINK companion availability, FAM structure,
BED dimensions, BIM identifier conventions, and gzip BED4 annotations. A
read-only check uses `validate_once(paths, contract, operation)` only when its
contract contains the validator identity and all relevant options. Reuse
requires unchanged file size, modification/change times, device, and inode;
the cache also checks identities after parsing, including resources inspected
by nested cached validators. Reused checks acquire each
consumer without rescanning. These are current-invocation identities, not
durable content fingerprints; checkpoint validation remains separate, and
sessions cannot be restored or reopened for another invocation.

Input-dependent checks—such as exact variant overlap between a newly formatted
table and a PLINK BIM—remain deferred until that table exists. Changed resources
are refused where recorded identity evidence is reused. The formatter similarly
reuses the validated VCF header and exact single-sample check, performs one full
extraction for the current VCF state, and creates every assigned module input
from that in-memory table. The method-specific compatibility rules remain with
their scientific consumers; the common layer does not change allele conventions,
sample-size semantics, missingness rules, or scientific thresholds.

`core.validation_reporting` groups observed checks into file sections using the
shared terminal styles; it does not read scientific inputs. Pipeline orchestration writes an
atomic YAML audit at `pipeline.validation.report_file`, first for startup and
again before announcing successful completion or when execution fails,
including observed later checks and failures. The
configured `logging.file_validation.max_screen_files` optionally limits passed
file sections (the default is 20), never problems or saved records. Additional
successful files are counted on screen and retained individually in the audit.
Direct commands use a record-only context with the same renderer and an atomic
audit, without activating pipeline caching or imposing new input requirements.
Startup deferrals remain explicitly labelled; a successful
pipeline is not used to invent passes for unrecorded checks. Unrecognised existing
files and validated inputs cannot be overwritten as validation reports.

This is an incremental consolidation, not a claim that every file now receives
full-content validation. Opaque resource checks and some existing resource
checks are availability-only; candidate-dependent checks and method-specific
stage summaries remain. See [Pipeline Input Validation](wiki/core/pipeline-input-validation.md)
for the current coverage and interpretation contract.

## Shared HTML table interaction

Every direct-module and pipeline HTML report is published through
`core.io.reports.write_html_report`. The writer adds one self-contained,
offline table controller: every column has a text filter and a sortable header,
column filters combine, and numeric values sort numerically with missing values
placed last. For JSON-backed paginated result tables, filtering and sorting use
the complete embedded result set rather than only the visible page. These
controls change presentation only; they do not modify, filter, or rewrite the
scientific result files.

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
