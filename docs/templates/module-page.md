# Module title (`postgwas command`)

> Documentation template: replace every instruction in angle brackets and
> remove this notice before publishing the page.

Use the opening summary to answer the reader's four questions in order:

> **What does this do? → What do I need? → How do I run it? → What do I get?**

## Purpose

<Explain in plain scientific language what the module produces and why a GWAS
researcher would use it.>

## What the analysis does

<Describe the implemented computational and scientific operations rather than
repeating function names. Explain scientifically material transformations and
assumptions. Cite primary publications and current official tool documentation.>

## When to use it

<Give practical use cases and state when the module is inappropriate or
unnecessary.>

## Input requirements

<List every user-supplied input for standalone execution. Distinguish artifacts
created by preceding pipeline steps from external resources.>

Document every applicable compatibility requirement:

- required columns or VCF fields;
- genome build and chromosome convention;
- population, ancestry, and LD-reference compatibility;
- variant identifiers and allele convention;
- sample-size definition and study design;
- file schema and external-software version.

Write `Not applicable` for a compatibility category that does not affect the
module. Explain the implemented failure behavior when required information is
missing.

## Command

State the exact public command name and whether standalone and pipeline modes
are implemented.

```console
postgwas command --help
```

## Direct mode

<Provide complete CLI-first commands for every distinct analysis and input
route. Include every required input/reference; show the expected table schema.
Group extended settings beneath the corresponding method. Label placeholders
clearly; do not use development-machine paths or omit arguments with ellipses.>

## Pipeline mode

<Provide complete commands for each supported pipeline route. Include the entry
VCF, external resources and required method choices; omit intermediates generated
by the planner. Explain the dependency chain and direct/pipeline input difference.
If standalone-only, say so explicitly instead of inventing a pipeline command.>

## Parameters

<!-- BEGIN GENERATED CONFIGURATION -->
<Generate defaults from canonical YAML and types or constraints from validated
schema metadata. Generate CLI option names from the current parser.>
<!-- END GENERATED CONFIGURATION -->

Use this table shape for information that cannot yet be generated:

| Parameter | Required? | Default | Description | Allowed values/notes |
|---|---:|---|---|---|

The Input requirements, Direct mode, Pipeline mode and Outputs headings are
stable navigation anchors used by the README. Keep one canonical module guide;
link specialised recipes instead of maintaining competing installation paths.

Do not copy defaults manually. Clearly distinguish required user inputs, CLI
overrides, user-configurable values, and internal algorithm invariants.

## Processing steps

<!-- BEGIN GENERATED PIPELINE -->
<Generate pipeline dependencies and order from the canonical registry and
planner.>
<!-- END GENERATED PIPELINE -->

<Describe the module's major internal operations in execution order. Explain the
rationale for scientifically important transformations.>

## Outputs

| Output | Description | Important columns or fields | Interpretation |
|---|---|---|---|
| <Exact implemented path or pattern> | <Purpose> | <Columns/fields> | <Meaning> |

Separate validated final artifacts from temporary or incomplete files. State
which metadata and resolved configuration records accompany the outputs.

## QC and logs

<Describe log locations, QC reports, warnings, skipped or filtered variants,
failure conditions, external-tool errors, and how successful completion is
verified.>

## Interpretation

<Explain how a researcher should interpret the main outputs and what conclusions
the analysis does not support on its own.>

## Common problems

| Problem or error | Likely cause | How to check | Solution |
|---|---|---|---|
| <Implemented symptom> | <Evidence-based cause> | <Relevant log or validation> | <Supported action> |

## Limitations

<State unsupported designs, incomplete validation, scaling limits, and known
scientific or software constraints.>

## Scientific references

<List primary method publications and current official tool documentation.
Record versions or access dates when upstream behavior may change.>
