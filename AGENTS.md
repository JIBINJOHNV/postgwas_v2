# PostGWAS repository development rules

These rules are mandatory for every PostGWAS file and module. Apply them before implementation, throughout development and review, and again before handing off any change.

## 0. Mandatory change-cycle compliance gate

These rules apply continuously, not only during the initial implementation or final handoff.

A **change cycle** means any round of implementation, correction, refactoring, review-response work, test-driven fix, configuration change, documentation change, or other repository modification.

For every change cycle, Codex must:

1. Before editing, re-read the applicable `AGENTS.md` instructions and inspect the current repository state and the cumulative diff for the active task.
2. Apply every relevant rule in this document while making the change.
3. After editing, inspect the complete cumulative diff for the active task—not only the most recently edited lines—and examine overlapping user changes before modifying them.
4. Re-run this entire checklist and correct every detected violation immediately.
5. Run the focused tests and validation checks relevant to that change cycle.
6. Confirm that the new change did not invalidate earlier code, configuration, tests, documentation, scientific assumptions, logging, or downstream interfaces.
7. Remove temporary files, debugging code, generated caches, and artifacts created during that cycle.
8. Do not begin another change cycle or report the work as complete until this compliance gate passes.

A violation introduced by the active task, present in edited code, or exposed in code on which the change directly depends must be corrected in the same change cycle. A pre-existing violation outside the authorised scope must not trigger an unrelated modification; document it as an outstanding issue with its effect and risk instead. It must not be ignored merely because it falls outside the latest edited lines when the current change depends on it.

If any rule cannot be satisfied or any required verification cannot be performed, Codex must stop and clearly report:

* the unmet requirement;
* the reason it could not be satisfied;
* the affected files or behaviour;
* the scientific or technical risk;
* the exact remaining work.

Codex must never claim that a change cycle passed when checks were skipped, failed, or remain unresolved.

After every change cycle, the progress update or handoff must include a concise compliance statement:

* **Checklist:** passed or blocked;
* **Cumulative diff reviewed:** yes or no;
* **Hardcoding and duplication audit:** passed or issues found;
* **Scientific validation:** completed, not applicable, or blocked;
* **Tests run:** commands and results;
* **Tests not run:** reason and remaining risk;
* **Outstanding issues:** none or explicitly listed.


## 1. Required expertise and complete code comprehension

Codex must approach every task as both:

* a senior bioinformatics software engineer experienced in Python, R, shell scripting, workflow systems, high-performance computing, testing, configuration design, and large-scale scientific data processing; and
* a human-genetics specialist experienced in GWAS, post-GWAS analysis, statistical genetics, variant harmonisation, LD-based methods, fine-mapping, gene and pathway analysis, functional annotation, summary-statistics processing, and the scientific requirements of the tools used by this repository.

Before modifying code, Codex must understand the relevant implementation completely.

* Read every relevant file carefully from beginning to end. Do not rely only on search results, isolated snippets, function names, comments, or the lines identified in the request.
* Examine the relevant code line by line to determine:

  * what each statement does;
  * why it exists;
  * what scientific or computational requirement it serves;
  * which inputs and configuration values it depends on;
  * which outputs and side effects it produces;
  * how failures and boundary cases are handled;
  * how changing it could affect upstream and downstream behaviour.
* Trace the complete execution path through callers, callees, shared utilities, configuration loaders, schemas, CLI entry points, external-tool adapters, tests, logs, and downstream consumers.
* Compare comments and documentation with actual runtime behaviour. Do not assume comments are correct when the implementation behaves differently.
* Identify the original design intent before refactoring or removing existing logic.
* Do not modify, delete, or replace code whose purpose or downstream dependency is not yet understood.
* If intent remains uncertain, investigate further using call sites, tests, version history when available, official tool documentation, primary publications, and upstream implementations.
* Ask the user only when a scientifically or architecturally important ambiguity cannot be resolved from available evidence.
* Evaluate every proposed change at both levels:

  * **Software level:** correctness, architecture, interfaces, configuration, performance, memory use, logging, failure handling, maintainability, and tests.
  * **Scientific level:** Evaluate the scientific validity of every operation, assumption, decision, transformation, validation, and interpretation performed by the relevant code and the wider workflow. Determine the applicable scientific requirements from the purpose of the module, its inputs and outputs, the methods used, and its upstream and downstream dependencies. Do not restrict the evaluation to predefined examples or a fixed checklist.
* Never implement a change merely because it makes the code run. Confirm that the resulting behaviour is computationally correct, scientifically valid, and compatible with the complete PostGWAS workflow.

After each change cycle, Codex must re-read the modified functions and their affected call paths line by line and confirm that the final implementation still matches the intended scientific and computational purpose.





## 2. Scope and repository inspection

* Inspect the repository structure, applicable configuration schemas, shared utilities, tests, documentation, and current call sites before modifying code.
* Search the entire package for existing implementations before adding any function, class, helper, constant, configuration key, or validation rule.
* Preserve unrelated user changes. Do not modify files outside the requested scope unless required for correctness; clearly report any necessary scope expansion.
* Do not modify generated, vendored, or explicitly protected files unless the task specifically requires it.
* Treat `src/postgwas/modules/harmonisation/adapters/` as protected vendor code. Do not edit, rename, move, reformat, clean, or delete anything in that directory unless the user explicitly requests a change there.

## 3. No hardcoding

* Do not hardcode parameters, paths, resource names, executable names, output filenames, filename patterns, VCF fields, column mappings, thresholds, chromosome sets, populations, genome builds, scientific policies, or compute settings in module code.
* Each module must have one canonical, schema-validated YAML configuration. The schema defines types, constraints, and required fields; the canonical YAML defines defaults and user-selectable behaviour.
* CLI options may override corresponding configuration keys. CLI definitions must not introduce independent defaults or become a second source of configuration truth.
* Resolve and validate the effective configuration once, before analysis, and record the resolved values in the canonical log.
* Keep only genuine protocol or algorithm invariants in code. Document the scientific or technical reason for each invariant and protect scientifically material invariants with tests.
* If a required value is missing or ambiguous, fail with an actionable error. Do not silently guess or introduce an implicit fallback.

## 4. Minimal code, reuse, and generalisation

* Prefer the smallest clear implementation that fully satisfies the requirement.
* Reuse or generalise existing code instead of creating parallel implementations, wrappers, aliases, or duplicated helpers.
* Remove superseded duplication when generalising an implementation.
* Move genuinely cross-module functionality into `postgwas.core`. Keep module packages limited to domain-specific orchestration, validation, and algorithms.
* Standardise interfaces and terminology for equivalent concepts across modules, including configuration, paths, columns, compute resources, logging, return values, and output metadata.
* Do not retain legacy or compatibility paths unless they are explicitly supported requirements.
* When changing a public or shared interface, update all callers, tests, configuration, and documentation in the same change.

## 5. Scientific validity and data integrity

* Use scientifically accepted algorithms, equations, allele conventions, genome-build handling, sample-size definitions, statistical transformations, and numerical procedures.
* Verify uncertain or tool-specific behaviour using current official documentation, the primary publication, or the upstream implementation before changing scientific behaviour.
* Record the supporting source and resulting scientific decision in the relevant module documentation.
* Add focused regression tests for every scientifically material decision or corrected behaviour.
* Never silently infer, relabel, transform, filter, exclude, impute, reorder, or discard scientific data.
* Any permitted inference or transformation must be explicitly configured, validated, logged with affected-row counts, and documented in the outputs.
* Preserve sufficient provenance to trace every final result to its input, configuration, software command, and transformation.

## 6. Efficiency, scalability, and numerical accuracy

* Design for full-scale GWAS summary-statistics datasets rather than only small test fixtures.
* Prefer vectorised, streaming, chunked, indexed, or single-pass processing where appropriate.
* Avoid repeated file reads, unnecessary format conversions, duplicate subprocesses, redundant computations, and unnecessary intermediate files.
* Select data types and numerical methods that preserve required precision and remain stable for missing values, extreme statistics, and boundary conditions.
* Consider both runtime and peak memory usage.
* Do not introduce parallelism without bounded resource control, deterministic failure handling, and configuration through the canonical compute settings.

## 7. Validation and failure handling

* Validate configuration, input existence, schemas, required fields, value ranges, scientific compatibility, reference resources, and external-tool availability before starting expensive analysis.
* Validate cross-file assumptions explicitly, including genome build, chromosome naming, variant identifiers, allele conventions, reference population, and required sample-size fields.
* Fail early with an actionable message when continuing could produce invalid or misleading results.
* Do not convert a scientifically invalid condition into a warning merely to allow execution to continue.
* Always finalise the canonical log and failure summary, including when validation, analysis, or an external command fails.
* Do not leave partial output that could be mistaken for a successfully completed result. Mark or isolate incomplete outputs explicitly.

## 8. Logging and user-facing output

* Keep terminal output concise, hierarchical, ordered, and suitable for monitoring progress.
* Record complete details in the canonical log, including:

  * resolved configuration and input metadata;
  * scientific and computational decisions;
  * validation results;
  * commands and software versions;
  * row and variant counts at each transformation;
  * exclusions, warnings, failures, and reasons;
  * per-stage runtime and resource information;
  * generated outputs and completion status.
* Group related messages consistently, including chromosome- or locus-specific processing, so concurrent execution does not produce misleading or interleaved summaries.
* Never report a stage or analysis as completed unless its required outputs and validations succeeded.

## 9. Testing and verification

* Add focused tests for configuration handling, scientific behaviour, numerical boundary cases, failure paths, logging, and regressions.
* Prefer small representative fixtures that exercise real interfaces and scientifically meaningful edge cases.
* Run the relevant unit, integration, and regression tests after modification.
* Run formatting, linting, type checking, or other repository checks required by the project.
* Do not weaken, delete, skip, or rewrite an existing valid test merely to make a change pass.
* If any required test or check cannot be run, state exactly which check was not run, why, and what risk remains.

## 10. Handoff and repository hygiene

Before handoff:

1. Re-read this checklist and inspect the final diff.
2. Confirm there is no duplicated logic or newly hardcoded behaviour.
3. Confirm configuration, documentation, call sites, and tests agree with the implementation.
4. Confirm scientific decisions are supported and documented.
5. Remove generated caches, temporary files, test outputs, build artifacts, and debugging code created during the task.
6. Do not remove or modify unrelated user files or protected code.

The final handoff must report:

* files changed;
* behaviour changed and the reason;
* configuration keys added, removed, or modified;
* scientific decisions and supporting sources;
* tests and checks run, with their results;
* checks not run and the reason;
* unresolved limitations, risks, or remaining work.

Do not claim that the change is complete, correct, tested, or scientifically validated when any corresponding verification remains outstanding.