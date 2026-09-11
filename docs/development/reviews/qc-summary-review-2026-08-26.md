# PostGWAS `qc_summary` — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/qc_summary/` (assessment.py 793, reporting.py 436, service.py 385, cli.py 111), `config/models/modules/qc_summary.py`, `config/defaults/modules/qc_summary.yaml`, `docs/wiki/modules/qc-summary.md`, and the module's pipeline wiring and contract with `harmonisation` and `filtering`.
**Method:** read-only on the tree — no file was modified. Behavioural claims about bcftools were verified against **bcftools 1.19 / htslib 1.19**; the performance claims in **I3** are measured, not estimated, on a synthetic 5,000,000-row table using this module's own column schema, null handling and expressions. Line numbers refer to the working tree as of 2026-08-26.

**Summary.** This is the strongest-engineered module reviewed so far: one VCF pass, atomic report publication, a completion manifest with a configuration digest that includes the bcftools version, and — unusually — accounting invariants that *raise* rather than warn. The QC arithmetic it performs is correct (see *Verified correct*). Two things undermine it. First, the module cannot tell a **field that does not exist** from a **value that is missing**, and with the shipped defaults that difference is the difference between a normal report and one that fails 100% of variants without erroring. Second, in pipeline mode it assesses the *filtered* VCF while every label, the JSON, and the log call it the raw merged VCF.

---

## CRITICAL

### C1 — An absent INFO/FORMAT tag is indistinguishable from a missing value, and the shipped defaults then fail every variant, silently

**Where:** `assessment.py:65` (`allow_undefined_tags=True`); `qc_summary.yaml:15, 18, 13` (the three `remove` actions), `:50` (`null_values: ['', '.']`); `assessment.py:270-298`.

**What.** The extraction deliberately runs `bcftools query --allow-undef-tags`. Verified against bcftools 1.19:

| input | rendered by bcftools | after `null_values: ['', '.']` |
|---|---|---|
| INFO tag not declared in the header | `.` | null |
| FORMAT tag not declared in the header | *empty string* | null |
| declared tag, no value for this variant | `.` / empty | null |

All three collapse to the same null. `_missing()` (`:122`) therefore cannot distinguish them, and every rule that consults a missing-value policy treats "this file has no frequency field" exactly as "this variant has no frequency".

The shipped defaults are `missing_af_action: remove`, `missing_info_action: remove`, `missing_pvalue_action: remove`. So if the VCF simply does not declare `INFO/<reference_af_column>` (default `EUR`), or `INFO/AF`, or `FORMAT/AF`, or `FORMAT/SI`:

```python
missing_for_comparison = missing_study_info_af | missing_external_af      # true on every row
"failure": discordant | (missing_for_comparison if af_missing_action == "remove" else …)
```

**every variant fails QC.** `qc_passed.num_records == 0`, `retained_fraction == 0.0`.

**Nothing raises.** The accounting invariants at `:570-581` balance perfectly — 0 passed + N excluded = N raw, and the rule matches reconcile — so `accounting_balanced` is written `True`, the completion manifest records `qc_passed_variants: 0`, and the log records `STATUS=COMPLETED`. The user receives a clean, internally consistent report stating that none of their variants passed quality control.

**Reachability.** Harmonisation does annotate `INFO/AFR,INFO/EAS,INFO/EUR,INFO/SAS` (`vcf_processing.py:549`), so a VCF straight out of this pipeline normally carries the default `EUR` tag. The exposed paths are: a study pointed at `--reference-af-column` for a population its reference file did not carry; any GWAS-VCF not produced by this pipeline (a documented standalone use); and a harmonisation run whose EAF-annotation step was skipped or produced an empty reference.

**Fix — the helper already exists.** `declared_vcf_tags(header, "INFO")` and `declared_vcf_tags(header, "FORMAT")` (`core/vcf.py:38`) return the IDs a header actually declares. Read the header once before extraction, check the ten configured fields against it, and report an undeclared field as its own named condition — fatal, or reported and excluded from rule evaluation — rather than folding it into per-variant missingness. Keep `--allow-undef-tags` afterwards: the point is not to make bcftools fail, it is to say *which* field is absent. This single check separates "3% of this study's variants have no frequency" from "this file has no frequency field at all", and those two currently produce the same report.

### C2 — In pipeline mode the module assesses the filtered VCF and labels every number "raw merged VCF"

**Where:** `planner.py:119-120`; `runners.py:171, 194, 335`; `cli_overrides.py:33`; `service.py:186-207`; labels at `reporting.py:109, 185`, `assessment.py:583-587`, `:786`.

**What.** Three facts compose:

1. The planner appends `qc_summary` **last** — after `sumstat_filter` (`:91`), `post_imputation_filter` (`:97`), `annot_ldblock` (`:100`) and `formatter` (`:111`).
2. Those runners mutate the shared namespace in place: `args.vcf = outputs["filtered_vcf"]` (`runners.py:194`), `= outputs["annotated_vcf"]` (`:171`), `= outputs["GRCh37"]` (`:335`).
3. `explicit_overrides` treats **any attribute present on the namespace** as an explicit CLI override — `if hasattr(namespace, destination)` (`cli_overrides.py:33`) — so the mutated `args.vcf` silently becomes `modules.qc_summary.inputs.vcf`.

`run_qc_summary_runner` (`runners.py:606-615`) is the only runner in the file that never consults `ctx`; it simply calls the direct entry point and inherits whatever the previous step left behind.

Every label then misstates what was measured: the screen report prints `2. Raw merged VCF` and "the raw merged VCF is unchanged"; the JSON `definition` says "every raw merged-VCF record"; the log records `vcf_output="raw_merged_only"`; and the docs say "Use QC summary after harmonisation".

**Scientific consequence.** `retained_fraction` is computed against a file that the same thresholds already produced, so it approaches 1.0 by construction. A reader sees *"QC-passed variants: 8,412,336 (99.94% of raw VCF)"* and concludes the study is clean — while the real losses happened in `sumstat_filter` and are invisible here. That defeats the module's stated purpose, which the docs give as "understand the effect of candidate filtering rules".

It also produces **two** QC assessments per pipeline run that claim to describe the same thing. Harmonisation calls the identical assessment on the genuine merge (`harmonisation/service.py:4671`, `vcf_path=raw_vcf_path`), so the user gets one true "raw merged VCF" report from harmonisation and a second, differently-numbered "raw merged VCF" report from the `qc_summary` step.

**Fix.** Decide what QC is contracted to describe and make it explicit instead of inherited. Either take the harmonised merge from `ctx` in `run_qc_summary_runner` (every other runner already reads `ctx`, and `get_vcf_from_context` exists at `runners.py:115`), or order QC before `sumstat_filter` in the planner. If assessing the processed file is genuinely intended, then stop calling it raw and record in the report which upstream steps produced the input.

---

## IMPORTANT

### I1 — `qc_summary` and `filtering` are two implementations of one QC policy, and their shipped defaults already disagree

**Where:** `config/defaults/modules/qc_summary.yaml:11-34` vs `config/defaults/modules/filtering.yaml:6-28`.

| Rule | `modules.qc_summary.rules` | `modules.filtering` |
|---|---|---|
| MAF floor | `maf_min: 0.01` | `maf_min: 0.01` |
| INFO floor | `info_min: 0.7` | `info_min: 0.7` |
| **INFO ceiling** | **`info_max: 1.05`** | **`info_max: null` — no upper bound** |
| missing AF / INFO / P action | `remove` | `remove` |
| AF-difference cutoff | `maximum_af_difference: 0.2` | `frequency_difference_max: 0.2` |
| reference tag | `reference_af_column: EUR` | `reference_population_tag: EUR` |
| palindromic | 0.4–0.6, remove | 0.4–0.6, remove |
| MHC | full GRCh37/GRCh38 interval table | the same table, written again |
| Neff outliers | `sample_size_outlier_standard_deviations: 5.0` | absent |

Same science, two schemas, two names for each of two identical quantities, the MHC coordinate table duplicated verbatim, and one default that **already diverges**. The docs say QC exists so users can "understand the effect of candidate filtering rules" — with the shipped defaults it does not, because a variant with INFO 1.2 fails the QC preview and survives the actual filter.

The codebase already knows this class of problem exists: `harmonisation/service.py:4660-4670` emits a purpose-built warning when `info.mach_rsq_max` disagrees with `modules.qc_summary.rules.info_max`. That guard was written for one pair of keys; the `qc_summary` ↔ `filtering` pair has none.

This is the largest reusable-function opportunity in the module: one shared rules model and one MHC interval table consumed by both, with `filtering` materialising a VCF and `qc_summary` only counting. Short of that, reconcile `info_max` deliberately and add the same cross-check warning.

### I2 — The effective-sample-size outlier check is one-sided, and it is the wrong side

**Where:** `assessment.py:385-387, 418-421`; `reporting.py:213-224, 379-387`.

The module computes only an upper bound and counts only above it:

```python
outlier_threshold = pl.col(mean_column) + sample_size_sd_multiplier * pl.col(sd_column)
_count(stage_sample_size & (sample_size > outlier_threshold), prefix + "__…above_outlier_threshold")
```

In GWAS summary statistics a **high** Neff is not a defect. A **low** Neff is: variants carried by only a subset of contributing cohorts, or by poorly imputed ones, carry inflated and unstable effect estimates, and are the standard target of sample-size QC (Winkler et al. 2014, *Nat Protoc*, recommend removing variants whose N falls below a fraction of the study maximum).

With the default `k = 5.0` and a realistically tight Neff distribution the upper count is structurally near zero. In my 5M-row simulation with Neff ~ N(150 000, 4 000) it returned **0** in both the raw and QC-passed stages. The module spends a mean, an SD, a threshold and two counts per stage to report a number that is almost always zero, and reports nothing about the tail that matters. The minimum Neff *is* reported (`:405-407`), but no rule, warning or threshold is attached to it.

Add the lower tail. Keep the upper count if it is wanted for symmetry, but the one that changes a decision is the low one.

### I3 — `.over(pl.lit(1))` defeats the streaming aggregation: measured 2.3× peak memory and 1.4× runtime

**Where:** `assessment.py:425-443`.

`_with_sample_size_statistics` broadcasts each stage's Neff mean and SD back onto every row with a window over a constant, adding four Float64 columns of full table height before the aggregation runs.

Measured on a 5,000,000-row / 316 MB table with this module's column schema, `null_values`, casts and mask logic (polars 1.44.1):

| implementation | wall time | peak RSS |
|---|---|---|
| current — `.over(pl.lit(1))` then aggregate | 6.3 s | **899 MB** |
| two-pass — aggregate the four scalars, then `pl.lit(...)` | 4.4 s | **396 MB** |

Results identical to float precision (`raw_mean` 149999.9057102168 vs 149999.90571022; `raw_sd` 4000.498526033483 vs 4000.49852603352). The gap scales with rows — 250 MB vs 157 MB at 1M, 899 MB vs 396 MB at 5M — roughly **+100 MB of avoidable peak per million variants**. A 20M-variant meta-analysis pays about 2 GB for it.

The docs' Processing step 5 claims "streaming Polars aggregation to calculate … in one pass". The window broadcast is what prevents that from being true.

**Fix.** A small first `select` of the four scalars, collected streaming, then `pl.lit(...)` in `_metric_expressions`. Two streaming scans of a temporary TSV cost far less than one materialised scan.

**Do not** take the tempting shortcut of deriving the SD from sum and sum-of-squares in a single pass. Neff is nearly constant within a study, which is precisely the case where `E[X²] − E[X]²` loses most of its significant digits to cancellation.

### I4 — `--threads` is accepted, mapped into the configuration, and never applied

**Where:** `service.py:50` (`"threads": "execution.threads"`); nothing in `modules/qc_summary/` reads it.

Polars sizes its thread pool at import time from the core count, so `postgwas qc --threads 4` on a 64-core shared node still spawns 64 worker threads. On a scheduler that enforces a CPU allocation this is the classic oversubscription failure.

This is the same defect the 2026-08-25 harmonisation review raised as its Critical #5, and **the fix already exists in this repository**: `harmonisation/service.py:221-224` sets `POLARS_MAX_THREADS` in the parent before any pool starts, with a comment recording that polars treats it as process-start configuration. `qc_summary` is a separate process entry point, so the same export applied at the top of `run_qc_summary_direct` — before polars is imported — is sufficient. Grep confirms `POLARS_MAX_THREADS` currently appears exactly once in `src/`.

### I5 — The VCF's genome build is never checked against `target_build`, and `remove_mhc` is on by default

**Where:** `service.py:208`; `assessment.py:334-361`; `qc_summary.yaml:8, 24`.

`genome_build = module.target_build.value` selects the MHC interval and names every output file. Nothing reads the VCF header. The shipped default is `target_build: GRCh37` with `remove_mhc: true`, and `--genome-build` is declared with `suppress_default=True` (`cli.py:44-48`), so omitting it silently accepts GRCh37.

Assess a GRCh38 VCF that way and the MHC rule excludes `chr6:28,477,797–33,448,354` interpreted as GRCh38 coordinates — an interval offset by roughly 32 kb from the true GRCh38 MHC. The "variants inside the configured MHC region" count is then wrong at both edges, and the QC-passed subset excludes the wrong variants, with no error.

`validate_vcf_header_contract` (`core/vcf.py:122`) already parses the declared build out of the `##genome_build=` header and is already used by `filtering` (`sumstat_filter.py:395`) and `ld_clumping` (`service.py:401`). `qc_summary` is the third VCF consumer that skips it. Call it with `required_fields=()` so the field-level reporting C1 describes stays in this module's hands, and compare `declared_build` to `target_build` exactly as `ld_clumping/service.py:414-419` does.

---

## OPTIONAL

**O1 — `INFO/AF` and `FORMAT/AF` are the same number, extracted and reported twice.** The GWAS-VCF writer sets both from one value: `record.info["AF"] = result.alt_freq` (`adapters/gwas2vcf/vcf.py:220`) and `record.samples[trait_id]["AF"] = result.alt_freq` (`:231`), where `alt_freq` is the study's own effect-allele frequency (`gwas.py:333-337`). Harmonisation's `population_frequency_qc.study_field: '%INFO/AF'` confirms INFO/AF is the study frequency, not a reference one, and nothing rewrites either tag after VCF creation. So `study_info_af` and `study_format_af` are duplicate columns in the ten-field extraction, `study_af_missing` and `format_af_missing` are duplicate metrics, and `missing_af_action` is applied twice to one underlying value. The module's own fixture writes them equal on every row (`tests/test_qc_summary_af.py:36-43`). If tolerating VCFs where they differ is intentional, document why; otherwise extract one.

**O2 — the section-renumbering hack in `qc_summary_lines`.** For the standalone path the function builds the whole "1. Before VCF creation" block, then finds it by searching rendered output for the literal strings `"1. Before VCF creation"` and `"2. Raw merged VCF"` (`reporting.py:405-413`), deletes the slice, and runs three `str.replace` calls across **every** line to renumber (`:415-426`). The `next(...)` calls raise a bare `StopIteration` if the rendering ever changes, and the blanket replace rewrites those substrings wherever they appear, including inside a filename or a rule criterion. Emitting the pre-VCF block only when `include_pre_vcf` and numbering from a counter is shorter, faster, and cannot mis-fire.

**O3 — a resumed run prints a different report than a fresh one.** `_write_assessment_reports` serialises the assessment to JSON *before* `assessment["reports"]` is assigned (`assessment.py:765`), so the on-disk JSON has no `reports` key. On resume, `_load_completed_assessment` returns that JSON and `qc_summary_lines` falls back to `reports = assessment.get("reports") or {}` (`reporting.py:98`), silently dropping the "Metric report" and "Rule report" lines a fresh run prints. Assign `reports` before serialising.

**O4 — `accounting_balanced` is self-certified.** It is written as the literal `True` (`assessment.py:594`) only after the checks pass, and `_load_completed_assessment` then gates resume on reading it back (`service.py:136`). It can only ever be `True` in a file this code wrote, so it adds nothing over the checksum manifest that already guards resume.

**O5 — Ts/Tv is reported with no expected range.** It is the most diagnostic sanity metric in the report — genome-wide GWAS should land near 2.0–2.1, and a value near 0.5 means the allele columns are scrambled — yet it renders as neutral "analysis" text (`reporting.py:196, 346`) while far less informative counts get warning colouring through `metric_warning_kind`. A configured plausible interval with a warning would cost a few lines and catch a whole class of upstream corruption.

**O6 — `_all_rows()`** (`assessment.py:147-150`) builds `CHROM.is_null() | CHROM.is_not_null()` to obtain an always-true mask. `pl.lit(True)` is the same value without touching a column, and the comment explaining the trick becomes unnecessary.

**O7 — the polars floor in `pyproject.toml:43` is wrong.** The declared constraint is `polars>=0.20`, but the module calls `collect(engine="streaming")` (`assessment.py:518`). Verified against polars 0.20.31: `LazyFrame.collect` has no `engine` parameter there, and the argument is **silently absorbed rather than rejected** — so an environment that resolves to 0.20 runs the whole aggregation non-streaming with no warning, which is exactly the failure mode I3 is trying to remove. Raise the floor to a version that implements `engine=`.

**O8 — every rule mask is built twice**, once for its `rule_NN__failed_raw` count and again inside the `combined_failure` OR chain (`assessment.py:478-492`). Polars' common-subexpression elimination may collapse them. Measure before changing anything — this is a note, not a defect.

---

## Verified correct — do not re-flag

- **MHC intervals.** GRCh37 `chr6:28,477,797–33,448,354` and GRCh38 `chr6:28,510,120–33,480,577` are the GRC region definitions for their builds, and the docs cite the GRC pages directly. The comparison strips a `chr` prefix from both the configured value and the data (`assessment.py:336-344`), so `chr6` and `6` both match, and the interval is applied inclusively as documented.
- **The MAF rule is correctly two-sided on an ALT-frequency column** — `af < maf | af > 1 − maf` (`:229`) — which is the right way to impose a MAF threshold on an effect-allele frequency.
- **The palindromic band [0.40, 0.60]** matches harmonisation's, already cleared in the 2026-08-25 and 2026-08-26 harmonisation reviews.
- **Ts/Tv construction.** Transversion is defined as `snp & ~transition` rather than enumerated (`:178`), so the two classes are exhaustive over SNPs and the ratio cannot be distorted by an unclassified remainder. Non-SNPs are excluded from both.
- **Per-variant missingness detection is sound.** `table.null_values: ['', '.']` covers both what bcftools emits for an absent INFO value (`.`) and for an absent FORMAT value (empty) — verified against bcftools 1.19. What it cannot do is separate that from an absent *field* (C1).
- **The accounting invariants are a genuine guard** (`:570-581`): raw = excluded + passed, and total rule matches = unique failures + overlapping matches, both raised as errors. Note only what they cannot catch — they hold perfectly when every variant fails.
- **Report publication is atomic.** Reports are staged into `NamedTemporaryFile` handles in the destination directory and `.replace()`d (`:617-697`); the temporary extraction table is removed in a `finally` whether the assessment succeeds or fails (`:789-793`).
- **The reuse claim in the docs is real, not aspirational.** `harmonisation/qc_reporting.py:211-222` is a thin wrapper over `qc_summary_lines`, and harmonisation calls the same `run_qc_assessment` (`harmonisation/service.py:161, 4671`). One policy, one implementation, exactly as documented — this is the model the `filtering` duplication in I1 should follow.
- **Provenance and resume machinery** — completion manifest, configuration digest that includes the resolved bcftools path *and* version, resolved-configuration export, structured `PARAM`/`INPUT`/`OUTPUT`/`STATUS` logging — is the strongest of the modules reviewed so far and needs no change.

---

## Documentation

`docs/wiki/modules/qc-summary.md` is thorough and mostly accurate. Four corrections:

- **Processing step 5** claims "streaming Polars aggregation … in one pass". Measured behaviour is a materialised window broadcast at 899 MB peak for 5M rows (I3).
- **"Use QC summary after harmonisation"** (When to use it) contradicts the planner, which runs it last, after filtering (C2). Whichever way C2 is resolved, these two must agree.
- **Common problems, row 1** — "Required VCF field is absent … Read preflight and bcftools diagnostics in the QC log" — describes a preflight that does not exist. `--allow-undef-tags` guarantees bcftools emits no diagnostic for an absent field, so there is nothing in the log to read (C1).
- Nothing states that **a 0% pass rate is the expected symptom of a missing reference-AF tag**. That is the most likely support question this module will generate, and one sentence in *Common problems* would answer it.

The *Interpretation* section correctly warns that concordance is meaningful only under matching orientation, build and population. Worth extending with the fact that the module verifies none of the three.

---

## Suggested order of work

1. **C1** — a header-declared-tag preflight built on `declared_vcf_tags`. Turns a silent, internally consistent 0%-pass report into a named error, and separates "field absent" from "value missing" everywhere downstream. ~15 lines.
2. **I5** — the build check via `validate_vcf_header_contract`, the same call two sibling modules already make; fold it into the same preflight as step 1.
3. **C2** — decide which VCF QC is contracted to describe, take it from `ctx` rather than from a mutated namespace, and make the labels tell the truth.
4. **I4** — set `POLARS_MAX_THREADS` before polars is imported, copying `harmonisation/service.py:221-224`.
5. **I3** — replace the window broadcast with a four-scalar first pass. Measured: −56% peak memory, −30% runtime, identical numbers.
6. **I2** — add the lower Neff tail.
7. **I1** — unify the rules model and the MHC table with `filtering`, and reconcile `info_max` deliberately rather than by accident.
8. **O1–O8.**

Steps 1–3 stop the module from reporting a confident wrong number. Steps 4–5 are pure cost with no scientific change. Steps 6–7 are the science this module should be doing and currently is not.
