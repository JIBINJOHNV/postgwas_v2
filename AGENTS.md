# PostGWAS repository development rules

These rules are mandatory for every PostGWAS file and module. Apply them before implementation, throughout development and review, and again before handing off any change.

## Global rule: missing required arguments

This rule applies to every public direct-module command and every pipeline command, in every supported execution mode.

* Resolve command-line overrides and canonical YAML first, then validate all mode-specific requirements before starting analysis, invoking external tools, or creating scientific outputs.
* If any required argument or resolved configuration value is missing, stop immediately with a non-zero exit status. Do not continue with a warning, guess a value, or use an undocumented fallback.
* The terminal and canonical log must use an actionable message that names the public option and, when applicable, its canonical configuration key: `Required argument not provided: --option-name. Provide --option-name VALUE or set configuration.path in the run configuration.`
* Report all independently missing required arguments together so the user can correct them in one attempt. Keep conditional requirements limited to the selected command, execution mode, analysis, and validated pipeline plan.
* Report a provided but nonexistent, unreadable, malformed, or scientifically incompatible value as an invalid argument, not as a missing argument.
* Required-argument detection and message formatting must use shared CLI/configuration validation rather than duplicated module-specific checks. Pipeline-generated intermediate inputs must not be required from the user.
* Regression tests must cover every registered public direct-module command and representative single- and multi-module pipelines. They must verify the non-zero exit, exact missing-option names, canonical-log entry, absence of analysis or external-tool execution, and correct handling of requirements satisfied through YAML.

## Global rule: CLI-help default labels

This rule applies to every displayed default on global, configuration, resource, validation, direct-module, and pipeline help pages.

* Render the label exactly as `Default:`. Do not qualify it as `Configured default:`, `Configured YAML default:`, `Effective default:`, `Packaged default:`, `Upstream default:`, or any other variant.
* Use the shared CLI-help formatter for the label and value. Module and pipeline parsers must not define, parameterise, or style their own default labels.
* The shorter label changes presentation only: the displayed value must still come from the canonical schema-validated YAML or the resolved effective configuration, and argparse must not become a second source of truth.
* Regression tests must inspect every registered public help page and every pipeline target, and must fail if any qualified, differently cased, or otherwise non-standard default label is displayed.

## Global rule: CLI-help option and description alignment

This rule applies to global, configuration, resource, validation, direct-module, and pipeline help pages at every supported terminal width.

* For every CLI option, display the option invocation and the beginning of its help description on the same row whenever both fit within the terminal width. Start every help description at the shared description column used by the surrounding options.
* If a help description wraps onto additional rows, indent every continuation row to that same shared description column. Never restart wrapped text beneath the option invocation.
* Place the entire help description on the next row only when the option invocation itself reaches or exceeds the shared description column. The description must still begin at the shared description column, and all continuation rows must retain that alignment.
* Treat `Required:`, explanatory text, `Default:`, available choices, conditions, and restrictions as parts of the help description and apply the same alignment rules to all of them.
* Implement this presentation through the shared CLI-help formatter. Module and pipeline parsers must not insert manual spaces, padding, or leading newlines to force alignment.
* The presentation rule must not change parsing, requiredness, resolved configuration, defaults, choices, usage text, or runtime behaviour.
* Regression tests must cover direct-module and context-sensitive pipeline help at representative narrow and wide terminal widths. They must verify same-row descriptions when space permits, aligned wrapped continuation rows, and aligned next-row descriptions for option invocations that do not fit before the shared description column.

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

## 8. Mandatory live progress reporting

* Every public scientific direct-module execution and every pipeline execution must show a shared live progress bar for the entire operation whenever terminal display is active. This top-level progress is mandatory and must not be disabled by `logging.show_progress`; that setting may control additional subordinate detail only. `--hide-screen` may intentionally suppress physical terminal display, but the durable progress summary must still be written to the canonical screen log.
* The progress display must remain visible and refresh while work is active, including while output is being captured by the shared screen recorder and while warnings or external-tool messages are emitted. It must show the active operation, completed operations, total operations, percentage, and elapsed time. Do not invent intermediate percentages when the underlying module exposes no measurable work units; keep the operation active and update elapsed time until its validated completion.
* Pipeline progress must follow the validated execution plan in order and update after each validated pipeline step. Modules with scientifically meaningful internal stages may show subordinate progress without changing top-level or pipeline completion accounting. Shared orchestration must provide the mandatory top-level progress for modules that do not expose internal stages; do not add duplicated module-specific progress implementations.
* A failed, interrupted, or non-zero operation must stop below 100%, name the failed operation or stage, and never be rendered as completed. A successful operation may reach 100% only after its required outputs and validations succeed.
* Progress presentation must use the shared `postgwas.core.ui` and screen-recording infrastructure, add negligible computational overhead, remain deterministic and safe around parallel workers and external tools, and be regression-tested for direct modules, pipelines, success, failure, terminal recording, hidden-screen logging, and nested detailed stages.

## 9. Logging and user-facing output

* Keep terminal output concise, hierarchical, ordered, and suitable for monitoring progress.
* Every public scientific direct-module CLI and every pipeline CLI help page must expose the centrally defined `--resume`/`--no-resume` and `--overwrite` controls. `run.resume: true` and `run.overwrite: false` are the canonical YAML defaults; CLI actions must use suppressed defaults so YAML remains authoritative, help must render both defaults with the shared green default style, and overwrite must take precedence over resume. Do not duplicate these options in module-specific parsers, and protect the complete registered command and pipeline surfaces with regression tests.
* Every public module and pipeline help page must end with clear, copy-ready command examples for each supported execution mode and scientifically distinct analysis offered by that command. An example must include every argument required for that exact mode, use the real public option names, show placeholders that make the expected file or value type obvious, and run after the user replaces only those placeholders. Direct-module examples must include their module-specific required inputs. Pipeline examples must include the pipeline entry inputs and every external reference or resource that the selected modules require but preceding pipeline steps do not generate; they must omit module-specific intermediate files supplied automatically by the pipeline. When requirements vary by analysis, input type, or configuration choice, provide a separate labelled example for each variant instead of an incomplete generic command. Include the canonical configuration-export command when the command supports run configuration. Keep example wrapping, indentation, section order, terminology, and styling consistent through shared CLI-help utilities rather than duplicating formatted help text in individual parsers.
* Every option whose value is mandatory for the displayed execution mode must begin its help description with the shared bold-red `Required:` label. Requiredness is mode-specific: direct help must mark every input that the module itself needs, while pipeline help must mark pipeline entry inputs and external resources that no preceding step generates, hide pipeline-generated intermediate inputs, and leave genuinely optional inputs unmarked. A mandatory value may be supplied through canonical YAML or an explicit CLI override unless the interface documents a scientifically necessary CLI-only exception; the visual marker must not create a second default or bypass resolved-configuration validation. Required-option metadata and styling must be centralised and protected by direct-module and context-sensitive pipeline regression tests.
* CLI-help regression tests must cover every registered public module and representative single- and multi-module pipeline help page. Tests must verify that displayed examples use accepted options, contain all mode-specific required arguments, exclude pipeline-generated intermediate inputs, remain consistent with canonical YAML and parser validation, and can be parsed after deterministic test placeholders are substituted. A help page containing an incomplete, obsolete, non-runnable, or misleading example is a release-blocking defect.
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

## 10. Global validated checkpoint, resume, and restart policy

* Apply the canonical global `run.resume_policy` to every scientific direct module and every pipeline. Do not create module-specific resume defaults or an independent policy outside the schema-validated global YAML.
* With default resume enabled, reuse a completed stage only after its checkpoint identity, resolved parameters, tracked input fingerprints, software identity, declared output paths, output fingerprints, and scientific output validation all match the current run.
* For a real partial run, resume at the next stage only when every earlier stage has a valid checkpoint and validated required outputs. Never infer completion from file existence, log text, timestamps, or a warning alone.
* When no valid sub-stage checkpoint exists, restart from the earliest safe module or pipeline-stage boundary. Never pass arbitrary partial scientific files to a later stage.
* When resolved parameters or tracked inputs changed and the global policy is `warn_and_restart`, emit a prominent terminal warning and canonical-log record, invalidate downstream checkpoints, and restart automatically. Delete or replace only outputs recorded as PostGWAS-owned, checksum-validated artifacts inside the configured output root; refuse automatic deletion for symlinks, paths outside that root, unknown artifacts, or externally modified outputs.
* `--overwrite` remains an explicit forced restart and takes precedence over resume. An automatic policy restart must be recorded as a restart decision; it must not silently mutate the resolved `run.overwrite` value.
* Treat `run.resume`, `run.overwrite`, `--resume`, `--no-resume`, and `--overwrite` only as execution controls: they must not contribute to the scientific checkpoint content digest or be restored from an earlier pipeline stage. A successful explicit overwrite must therefore be reusable by the next default resume. All content-determining resolved parameters remain checkpointed and must trigger the configured warning-and-restart policy when changed.
* Write each checkpoint atomically only after the stage outputs and scientific invariants succeed. Record the resolved configuration digest, input and output fingerprints, software versions, stage identity, metrics, and upstream checkpoint dependencies needed to prove that reuse is safe.
* Keep checkpoint manifest paths and behavior configuration-driven, standardise the implementation in `postgwas.core`, and add regression coverage for completed resume, partial-stage resume, changed-parameter restart, changed-input restart, externally modified output refusal, downstream invalidation, direct-module execution, and pipeline execution.

## 11. Testing and verification

* Add focused tests for configuration handling, scientific behaviour, numerical boundary cases, failure paths, logging, and regressions.
* Prefer small representative fixtures that exercise real interfaces and scientifically meaningful edge cases.
* Run the relevant unit, integration, and regression tests after modification.
* Run formatting, linting, type checking, or other repository checks required by the project.
* Do not weaken, delete, skip, or rewrite an existing valid test merely to make a change pass.
* If any required test or check cannot be run, state exactly which check was not run, why, and what risk remains.

## 12. Handoff and repository hygiene

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
