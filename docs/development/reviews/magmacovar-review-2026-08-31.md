# PostGWAS `magmacovar` — module review

**Date:** 2026-08-31 · **Scope:** `src/postgwas/modules/magmacovar/` (`service.py` 494, `main.py` 457, `cli.py` 206, `errors.py` 8), `config/models/modules/magmacovar.py`, `config/defaults/modules/magmacovar.yaml`, `docs/wiki/modules/magmacovar.md`, `tests/test_magmacovar.py`, and the contract with `magma` and the downstream consumer.
**Method:** read-only on the tree — no file was modified. Every constructed flag was checked against MAGMA's gene-property command line and every constrained vocabulary against MAGMA's own. Line numbers refer to the working tree as of 2026-08-31.

**Summary.** This is a small, tightly validated wrapper and it gets the precise things right: `direction-covar=` rather than the blanket `direction=`, MAGMA's own vocabulary in the enums (`smaller`, not `less`), missing-value handling delegated to MAGMA rather than reimplemented, and — the best thing here — a **pre-check of per-property missingness against MAGMA's own `max-miss` rule**, so a property MAGMA would silently drop is named and refused before the run. There are no Critical findings.

Two findings, both about consistency with the module's own siblings rather than about the MAGMA call: it applies **no multiple-testing correction** while `magma` ships a full correction framework operating on the identical file format, and the **two halves of a use case can only be selected separately** even though `single_cell` can select them as a unit.

**Implementation status (2026-08-31):** I1 is resolved. MAGMAcovar now writes a
separate, schema-validated corrected TSV with global Bonferroni, Šidák, Holm,
and BH-FDR values and a configured primary Bonferroni decision. The native
`.gsa.out` is preserved byte-for-byte as the FLAMES handoff. Both result files
and the native MAGMA log are completion-manifest artifacts. I2 remains open.

---

## IMPORTANT

### I1 — No multiple-testing correction, in the module built for large property screens

**Where:** nothing matches `bonferroni|fdr|adjust_p|multiple_testing|holm|sidak` anywhere in `src/postgwas/modules/magmacovar/`, `magmacovar.yaml`, or its config model. The output layout ends at `results_file: "{dataset_id}.gsa.out"`.

Its sibling does the opposite. `magma` ships a `multiple_testing` block with four global methods and four named gene-set families, and writes `*_gene_sets_corrected.tsv`. The two modules consume **the same file shape**: `magma`'s `correct_gene_set_p_values` reads `VARIABLE` and `P` from a `.gsa.out`; `read_magma_covariate_results` (`main.py:301-380`) validates exactly `VARIABLE, TYPE, NGENES, BETA, BETA_STD, SE, P` from a `.gsa.out` and already rejects duplicate `VARIABLE` values — which is precisely the precondition a correction needs.

**The scale is the point.** This module's own configuration cites Duncan et al. (2025) screening **461 cell-type properties** (`magmacovar.yaml:27`), and the shipped example runs GTEx tissue panels. A 461-test screen reported as raw p-values, when the sibling module would have Bonferroni- and BH-corrected the same table automatically, is a difference a reader will not expect between two commands from the same pipeline.

**It also flows downstream.** `ctx["magma_covar"]` becomes `args.magma_covariate_results_file` (`pipeline/runners.py:596`), and FLAMES carries a `magma_covariate_column` (`config/models/modules/flames.py:30`). Whatever that consumer does with the column, it receives uncorrected p-values from a screen of arbitrary width.

The docs do disclose it — "hypothesis family, and multiple-testing correction before interpretation" (`magmacovar.md:196`) — so this is a deliberate position, not an oversight. The objection is that the same position is not taken in `magma`, and a user running both gets corrected columns from one and raw p-values from the other for the same class of competitive test. Either correct here using the existing shared code, or say in both modules why the analyst owns it in one and not the other.

### I2 — The model and the direction of a use case can only be chosen separately

**Where:** `magmacovar.yaml:3-14` (`model_use_cases`), `cli/common.py:805-824` (rendered into help), `config/models/modules/single_cell.py:101` (`magmacovar_use_case`).

The configuration pairs a model with a direction deliberately:

```yaml
tissue_specificity:
  label: "Tissue specificity / FUMA-style"
  model: ["condition-hide=Average"]
  direction: greater
```

The pairing is not cosmetic — FUMA's tissue-specificity test is one-sided because only a positive association is interpretable there. `condition-hide=Average` with the shipped default `direction: two-sided` is a different analysis from the one the cited paper ran.

`model_use_cases` is consumed only to **render CLI help**: a "How to choose a MAGMA model" block that prints each use case's label, question, the `--covariate-model` string and the matching `--covariate-direction`, for the user to copy as two separate flags. Meanwhile `single_cell` names a use case in configuration (`magmacovar_use_case`) and gets both halves applied together.

So a magmacovar user who copies `--covariate-model condition-hide=Average` and does not also copy `--covariate-direction greater` silently runs a two-sided tissue-specificity test, and nothing warns. The mechanism to prevent that already exists one module over. Exposing the same selector — `--covariate-use-case tissue_specificity`, applying both fields, with the individual flags still available as overrides — is a small change that removes a real footgun.

---

## OPTIONAL

**O1 — one under-powered property fails the whole run.** `read_magma_covariate_results` raises when *any* COVAR row reports `NGENES < minimum_genes` (default 10). That is defensible — a property tested on fewer than ten genes is uninterpretable, and silently retaining it would pollute any correction — but on a 461-property screen a single sparse cell type blocks every other result. As property counts grow, excluding-and-reporting is likely the better contract than failing.

**O2 — a development note left in production code.** `# FIX: key is "magma"` at `pipeline/runners.py:450`.

**O3 — `--gene-covar` is never given a `use=` modifier.** MAGMA can restrict which covariate columns it reads at the file level; this module always reads the whole table and relies on `--model analyse=list,…` to choose what is tested. The result is equivalent, but a user with a 461-column table should know from the docs that every column is parsed and validated on every run.

---

## Verified correct — do not re-flag

**`direction-covar=`, not `direction=`.** MAGMA's `--model` offers `direction=`, `direction-sets=` and `direction-covar=`. For a `--gene-covar` analysis the covariate-specific modifier is the precise one, and it is what the module emits (`main.py:296`). The blanket form would have worked and been less exact.

**The direction vocabulary is MAGMA's, not English.** `MagmaCovarDirection = Literal["two-sided", "greater", "smaller"]` — `smaller`, which is what MAGMA accepts, rather than `less`. A test explicitly rejects a direction word MAGMA does not take (`test_magmacovar.py:140`).

**The command shape follows MAGMA's convention** (`main.py:262-297`):

```
magma --gene-results <genes.raw>
      --gene-covar <file> missing-values=<v> max-miss=<f> [missing-genes=fill]
      --model <modifiers…> direction-covar=<direction>
      --out <prefix>
```

with modifiers appended to the file argument as MAGMA expects. `--model` entries are rejected if they carry a leading dash, so a user cannot smuggle a second flag through the model list. The three `--gene-covar` modifier names are tied to the MAGMA manual in the docs (`magmacovar.md:42`) and pinned by tests (`:117`, and the fake-MAGMA run asserting the rendered command at `:217-220`); I did not have the MAGMA binary here to verify them independently, and the project's own citation plus test is the right way to hold them.

**Missing-value handling is delegated, not reimplemented.** PostGWAS passes `missing-values=median` and lets MAGMA impute. No second implementation of the statistics to drift.

**Per-property missingness is pre-checked against MAGMA's own rule** (`main.py:190-206`). MAGMA silently drops a property whose missingness exceeds `max-miss`; this module computes the fraction per property first, refuses, and names the property. A silently dropped tissue is a test that never happened and that the user would never notice — this is the module's most valuable check.

**The missingness denominator follows MAGMA's semantics exactly** (`:186-190`): `len(eligible_gene_ids)` when `missing-genes=fill`, because absent genes are then filled and count as missing, and `overlap` otherwise. Getting this backwards would either hide a failing property or reject a sound one.

**Constant properties are rejected** — fewer than two distinct finite values among overlapping genes — using a set that stops growing at two values, so the check costs nothing on a wide table.

**Input validation is thorough and located before any work**: `NA` as the only missing token (MAGMA's convention), duplicate gene IDs, non-numeric and non-finite values, ragged rows, duplicate header names, empty gene IDs, and non-UTF-8 files, each a named error with a line number. `maximum_missing_fraction` is bounded `le=0.2`, the ceiling the docs cite.

**Output parsing is defensive**: it filters to `TYPE == "COVAR"` — a `.gsa.out` carries other row types — requires the documented header, validates every statistic finite with `SE ≥ 0` and `P ∈ [0, 1]`, and rejects duplicate `VARIABLE` values.

**The upstream contract is enforced, not assumed.** `dependencies=("magma",)`, and `run_magmacovar_runner` obtains the gene results through `_validated_magma_gene_result`, which refuses unless the primary MAGMA mapping's `result_statistic_type` is `calibrated_gene_p_value`. An eMAGMA or chromMAGMA result that is not a calibrated gene p-value therefore cannot silently become the outcome variable of a gene-property regression.

**Staging, resume and ownership are handled**, with tests covering a zero exit that produced no output, resume revalidation without re-running MAGMA, overwrite preserving unrecorded files that merely share the output prefix, and refusing to remove unowned staging files.

**`model_modifiers` is the best in-config documentation in this codebase.** Each MAGMA `--model` modifier carries its syntax variants, whether it requires a value, a plain-language description, a worked example, **and a published citation** for how the modifier has actually been used. That is documentation that will still be true when the code changes.

---

## Documentation

`docs/wiki/modules/magmacovar.md` is accurate and unusually good at explaining *why* each model choice exists. Two additions:

- State plainly that **no multiple-testing correction is applied** and that the p-values in `{dataset_id}.gsa.out` are raw — particularly for wide screens, and particularly because the sibling `magma` module does correct (I1).
- State that a use case's **model and direction must be set together**, and that using the tissue-specificity model with the default two-sided direction is a different test from the published one (I2).

---

## Suggested order of work

1. **I2** — expose the use-case selector the `single_cell` module already uses, so model and direction cannot separate. Small change, removes a silent misconfiguration.
2. **I1** — decide the position on correction and make it consistent with `magma`. If corrections belong here, `correct_gene_set_p_values` already reads this exact file format.
3. **Documentation** — the two additions above.
4. **O1–O3.**

Nothing here changes a number MAGMA computes. I2 changes which test gets run without the user intending it.
