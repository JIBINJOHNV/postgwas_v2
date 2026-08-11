# Harmonisation — Correction Plan

Companion to `HARMONISATION_REVIEW.md`. This document says **what to change**, file by file.
Nothing here has been applied; the codebase is untouched.

Every "current code" block below was copied from the actual file and verified.

---

## How to work through this

| Phase | What | Effort | Why now |
|---|---|---|---|
| **P0-A** | Config defaults only — no code | ~30 min | Fixes 7 findings by editing one YAML file |
| **P0-B** | 12 code corrections | 1–2 days | These are the ones that produce wrong science or silent data loss |
| **P1** | 12 correctness/robustness items | ~1 week | Detection, provenance, resource validation |
| **P2** | Performance, structure, cleanup | ~1 week | Memory, resume, dead code |
| **P3** | Tests | ~3 days | Lock in P0/P1 so they cannot regress |

Do **P0-A first** — it is the cheapest risk reduction available and requires no code review.

---

# P0-A · Config defaults (single file)

File: `src/postgwas/config/defaults/modules/harmonisation.yaml`

Seven of these are cases where the file's **own `recommendation:` field already says the right
value** and the shipped `default:` says something else. Change the defaults to match.

| Policy | Current default | Change to | Why |
|---|---|---|---|
| `validation.max_reject_fraction` | `1.0` | `0.95` | At 1.0 the guard `fraction > fail_at` is arithmetically unreachable — `fraction` cannot exceed 1.0. The file's own `recommendation` is `0.95`. |
| `vcf.liftover_swap` | `exclude` | `update_tags` | `exclude` silently discards swapped variants **and** they are not counted as liftover failures (the YAML help says so). Own recommendation: `update_tags`. |
| `chromosome.rename_map` | `{23:X, 24:Y, 25:MT, M:MT}` | `{23:X, 24:Y, 25:XY, 26:MT, M:MT}` | In PLINK coding 25 is XY (PAR) and 26 is MT. Today 25→MT then gets dropped by `drop_mt`, and 26 is rejected. Own recommendation already states this. |
| `chromosome.allowed_after_split` | (differs from `allowed`) | make identical to `chromosome.allowed`, and add `XY` to both | MT is in `allowed` but not `allowed_after_split`, so MT is always dropped regardless of `drop_mt`. |
| `strand.af_discordance_action` | `warn` | `reject` | This is the only per-variant AF cross-check. At `warn` a transposed effect/other allele column produces a fully "successful" VCF with every effect direction inverted. |
| `population_frequency_qc.on_error` | `warn` | `fail` | The ancestry-mismatch detector currently fails silently and the run still reports OK. |
| `validation.declaration_mismatch_action` | `warn` | `fail` (at least for `effect.type`) | A wrong `odds_ratio` declaration makes step 03 compute `ln(ln(OR))` and reject ~50% of variants. |
| `population_frequency_qc.population_fields` | `{AFR, EAS, EUR, SAS}` | add `AMR: '%INFO/AMR'` | There are **five** 1000G super-populations. An admixed-American study can never match its own. |
| `vcf_processing.external_frequency_columns` | `[…, INFO/AFR, INFO/EAS, INFO/EUR, INFO/SAS]` | add `INFO/AMR` | Pairs with the above. |
| Palindromic strand policy | frequency distance | **Implemented:** full-study non-palindromic consensus; reference AF is post-orientation QC only | Between-population AF noise can no longer choose strand. |
| `position.min_value` | `0`, `recommendation: null` | `1`, and set `recommendation: 1` | Position 0 is the standard "could not be mapped" marker and is currently retained. The help text already says 1 is recommended but the structured field is null so nothing surfaces it. |

**Also correct two help texts that describe behaviour the code does not implement:**

- `input.null_values` claims *"Setting this key at all also replaces the per-chromosome read's own
  list at chromosome step 01."* It does not — `_chromosome_read_options` reads the separate
  `input.chromosome_null_values`. Either implement the inheritance or fix the sentence.
- `input.check_truncation` claims *"false — this preliminary stream check is skipped."* It is not
  skipped; `_inspect_text_file` runs unconditionally and only the *warnings* are gated.

---

# P0-B · Code corrections

## 1. Adapter exits with status 0 on every parameter error

**Severity: Critical** · Silent total data loss reported as success.

**Where:** `src/postgwas/modules/harmonisation/adapters/gwas2vcf/main.py`
lines 107, 110, 118, 125, 132, 142, 146, 150, 154, 158, 182

**Current:**
```python
        sys.exit()
```
Verified: `python3 -c "import sys; sys.exit()"` → **exit status 0**.

**Change to:**
```python
        sys.exit(1)
```
at all eleven sites. Line 182 is the column-index bounds check — the one guard that catches a
column-mapping off-by-one — so it matters most.

**And** in `src/postgwas/modules/harmonisation/gwas2vcf_runner.py:78`:

**Current:**
```python
    run_checked_command(
        cmd,
        "GWAS-to-VCF conversion for chromosome %s" % chromosome,
        logger=logger,
        error_type=RuntimeError,
        stdout_path=transcript,
        stderr_to_stdout=True,
        stdout_header="tool=gwas2vcf\ncommand=%s\n\n" % command_str,
    )
```

**Change to:** add the kwarg that already exists in `core/processes.py:130` and is used at
`vcf_processing.py:789`:
```python
    run_checked_command(
        cmd,
        "GWAS-to-VCF conversion for chromosome %s" % chromosome,
        logger=logger,
        error_type=RuntimeError,
        stdout_path=transcript,
        stderr_to_stdout=True,
        stdout_header="tool=gwas2vcf\ncommand=%s\n\n" % command_str,
        expected_outputs=["%s.gz" % output_vcf],   # <-- add
    )
```

**Also** `gwas2vcf_runner.py:89`:

**Current:**
```python
    return command_str, 0
```
The `0` is a literal, not the process status — and it is written into the QC JSON and the run
manifest as `gwas2vcf_exit_code`.

**Change to:** return the real status (requires `run_checked_command` to surface it), **or**
delete the field and its QC key rather than publish a fabricated success value.

---

## 2. `ALLELE1` / `A1` are mapped unambiguously to the effect allele

**Severity: Critical** · Inverts every effect direction for SAIGE-format inputs.

**Where:** `src/postgwas/config/defaults/modules/harmonisation.yaml`,
`module.sample_sheet_generator.column_aliases`, and
`src/postgwas/modules/harmonisation/sample_sheet_generator.py:161` (`_one_match`)

**Current (verified):**
```yaml
effect_allele_column: ['A1', 'AL1', 'ALLELE1', 'ALT', 'ALT_ALLELE', ..., 'EA', 'EFFECT_ALLELE', ...]
other_allele_column:  ['A2', 'AL2', 'ALLELE0', 'ALLELE2', 'OA', 'NEA', ..., 'REF', ...]
```

`ALLELE1` → effect is correct for BOLT-LMM and REGENIE (`ALLELE0`/`ALLELE1`) but **inverted for
SAIGE**, where `Allele1` is REF and BETA is aligned to `Allele2`. `A1`/`A2` means effect/other in
PLINK but first/second allele in several METAL-derived consortium formats.

`_one_match` raises only when two aliases *within one category* are present, so a file with
`Allele1`/`Allele2` produces a silent mapping — the repo's own test asserts
`result.warnings == ()` for exactly this (`tests/test_harmonisation_sample_sheet_generator.py:146`).

**Change to:**

1. Add a third alias tier alongside the ones that already exist for effect and frequency columns:
```yaml
ambiguous_allele_column: ['A1', 'A2', 'AL1', 'AL2', 'ALLELE1', 'ALLELE2']
```
2. Remove `A1`, `AL1`, `ALLELE1` from `effect_allele_column` and `A2`, `AL2`, `ALLELE2` from
   `other_allele_column`. **Keep `ALLELE0`** in `other_allele_column` — the `ALLELE0`/`ALLELE1`
   pair is unambiguous (REGENIE/BOLT), so match that pair as a *pair*, not as single tokens.
3. In `_row_for_file`, when an ambiguous allele token is matched: emit a mandatory review warning
   naming the two candidate conventions (PLINK/BOLT vs SAIGE) and set `requires_completion=True`
   (currently only set on missing EAF/sample size, `:568`) so the draft sample sheet cannot be
   consumed unreviewed.

---

## 3. MAF column adopted as EAF without the deferred swap flip

**Severity: Critical (config-gated)** · EAF describes the wrong allele for ~half of variants.

**Where:** `src/postgwas/modules/harmonisation/allele_frequency.py:1428-1441` and `:1588-1592`

**Current (verified):** the `inconclusive` + `accept` branch sets a provenance string and stops —
it never applies the `1−p` flip that the `else` (confirmed-`eaf`) branch applies at `:1443-1452`:
```python
                        if reference_decision == "inconclusive":
                            qc_info["maf_check_inconclusive_action"] = inconclusive_action
                            if inconclusive_action == "fail":
                                raise AlleleFrequencyError(...)
                            emit.warn(...)
                            qc_info["decision_source"] = "internal_eaf_reference_inconclusive"
                        else:
                            if STRAND_ACTION_COLUMN in cleaned_df.columns:
                                swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                                    "forward_swapped", "reverse_complement_swapped"
                                ])
                                cleaned_df = cleaned_df.with_columns(
                                    pl.when(swapped)
                                    .then(1.0 - pl.col(internal_eaf))
                                    .otherwise(pl.col(internal_eaf))
                                    .alias(internal_eaf)
                                )
                            qc_info["decision_source"] = "internal_eaf_reference_validation"
```

This matters because `strand.py:467` deliberately skips the EAF flip for a MAF-like column
(`… and not eaf_is_maf`) while `strand.py:452-453` has **already** rewritten EA/OA to the reference
orientation. The flip is *deferred*, not skipped — and on this branch it never arrives.

Compounding, at `:1588`:
```python
            frequency_is_aligned = (
                not study_says_maf
                or qc_info.get("maf_reference_decision") == "eaf"
            )
            if REFERENCE_AF_COLUMN in final_df.columns and frequency_is_aligned:
```
evaluates to `False` in exactly this branch, so the `strand_af_difference` concordance check — the
one gate that would notice — is skipped too.

**Change to:** hoist the flip out of the `else` so both branches apply it, and mark the column
aligned once it has been adopted:
```python
                        if reference_decision == "inconclusive":
                            qc_info["maf_check_inconclusive_action"] = inconclusive_action
                            if inconclusive_action == "fail":
                                raise AlleleFrequencyError(...)
                            emit.warn(...)
                            qc_info["decision_source"] = "internal_eaf_reference_inconclusive"
                        else:
                            qc_info["decision_source"] = "internal_eaf_reference_validation"
                        # applies to BOTH branches: the column is now being used as an EAF,
                        # so the swap flip strand.py deferred must be settled here.
                        if STRAND_ACTION_COLUMN in cleaned_df.columns:
                            swapped = pl.col(STRAND_ACTION_COLUMN).is_in([
                                "forward_swapped", "reverse_complement_swapped"
                            ])
                            cleaned_df = cleaned_df.with_columns(
                                pl.when(swapped)
                                .then(1.0 - pl.col(internal_eaf))
                                .otherwise(pl.col(internal_eaf))
                                .alias(internal_eaf)
                            )
```
and at `:1588`, treat "adopted as EAF" rather than "confirmed as EAF" as the condition:
```python
            frequency_is_aligned = (
                not study_says_maf
                or qc_info.get("maf_reference_decision") in ("eaf", "inconclusive")
            )
```

**Alternative, arguably better:** redefine `accept` to mean *"fall through to the external EAF
panel"* rather than *"use the unconfirmed column as-is"*. Then no flip is needed because the column
is not adopted at all.

---

## 4. Reference EAF check returns a positive answer when it has no power

**Severity: High** · A genuine MAF column is confirmed as an EAF and the confirmation is recorded
dataset-wide as `chromosome_reference_consensus`.

**Where:** `src/postgwas/modules/harmonisation/allele_frequency.py:705-710`
(`confirm_eaf_with_reference`)

**Current (verified):**
```python
    reference_minor_fraction = comparison["minor_fraction"]
    stats["reference_effect_allele_minor_fraction"] = float(reference_minor_fraction)
    if reference_minor_fraction > maf_cutoff:
        decision = "eaf"
    else:
        decision = "maf" if errors["maf"] < errors["eaf"] else "eaf"
```

If >95% of matched **panel** effect-allele frequencies are themselves ≤0.5, then panel-EAF ≈
panel-MAF and the two errors are identical — that is the *no-information* case, not evidence of
"eaf". For a modern imputed GWAS the matched set is dominated by rare variants, so this fires
routinely.

**Change to:**
```python
    reference_minor_fraction = comparison["minor_fraction"]
    stats["reference_effect_allele_minor_fraction"] = float(reference_minor_fraction)
    margin = float(pol.get("eaf.maf_decision_margin"))          # new policy, e.g. 0.02
    if reference_minor_fraction > reference_minor_cutoff:        # new policy, see below
        # The panel's own ALT is nearly always the minor allele, so comparing the study
        # column against panel-EAF and panel-MAF cannot distinguish them.
        decision = "inconclusive"
    elif abs(errors["maf"] - errors["eaf"]) < margin:
        decision = "inconclusive"
    elif errors["maf"] < errors["eaf"]:
        decision = "maf"
    else:
        decision = "eaf"
```

**Add two policies** (the cutoff is currently overloaded — see below):
- `eaf.reference_minor_fraction_cutoff`, default `0.95`
- `eaf.maf_decision_margin`, default `0.02`

**Related — `eaf.maf_decision_cutoff` is used for two unrelated quantities.** It is documented as
"the share of **study** frequency values at or below 0.5" (used at `:593` and
`validate_allele_frequency:535`), and then reused at `:707` as the cutoff on
`reference_effect_allele_minor_fraction` — a property of the **panel**. Tuning the study
MAF-detection sensitivity silently retunes the panel power gate in the opposite direction. The new
`eaf.reference_minor_fraction_cutoff` above resolves this.

---

## 5. Palindromic SNP resolution: implemented consensus-first policy

**Status: resolved.** PostGWAS no longer ranks the two palindromic orientations by study/reference
AF distance, so neither the old `&` ambiguity-band defect nor a winner-versus-runner-up margin can
select strand. The complete-study, non-palindromic evidence must meet
`strand.min_informative_variants` and `strand.consensus_threshold`; only that strong forward or
reverse consensus filters the palindromic candidate pool. Mixed, unresolved, or insufficient
consensus rejects the row as `palindromic_ambiguous`.

The aligned reference AF remains attached only for the post-orientation QC check. A difference above
`strand.af_tolerance` rejects a palindromic row by default through
`strand.palindromic_af_discordance_action`, but never changes `strand_action`, BETA, Z, or EAF.
Focused regression tests cover weak consensus, near-0.5 AF under strong consensus, and post-consensus
AF discordance.

---

## 6. The catastrophic-loss guard cannot fire

**Severity: High** · A study that loses 95% of its variants to a mis-detected scale completes with
status `OK`.

**Where:** `src/postgwas/modules/harmonisation/effect_validation.py:~496`

**Current (verified):**
```python
        fraction = (float(removed_total) / float(rows_in)) if rows_in else 0.0
        ...
        fail_at = pol.get("validation.max_reject_fraction")

        if rows_in and fraction > fail_at:
            raise EffectValidationError(...)
```
With `fail_at = 1.0` and `fraction ≤ 1.0` by construction, the branch is unreachable.

**Change to:** ship `0.95` (P0-A) **and** widen the denominator. The fraction is currently computed
from step-10 removals only; steps 03, 06, 07 and 08 each remove variants beforehand and none of
them count. `rows_read` is already available in `process_one_chromosome` at `:1253` — thread it
through and evaluate cumulative loss:
```python
        fraction = (float(rows_read - rows_out) / float(rows_read)) if rows_read else 0.0
```
Use `>=` rather than `>` so a configured `1.0` still catches total loss.

---

## 7. Real p-values below 1e-300 are deleted

**Severity: High** · Removes the lead SNPs of the largest GWAS whenever the study supplies BETA + P
and no SE — a very common public format.

**Where:** `p_values.py:717` → `standard_error.py:496-507` → `effect_validation.py:400-409`

**The chain:** any `0 ≤ p < pvalue.clip_low` (default `1e-300`) is flagged `__pval_clipped_low`;
that flag folds into `__se_from_clipped_pval`; `validation.se_from_clipped_pval` defaults to
`reject`.

But double precision holds p down to ~`5e-324`, so everything in `[5e-324, 1e-300)` is an exactly
representable, scientifically real p-value — **not** an underflow.

**Change (pick one):**

- **Preferred** — separate the two conditions. Only *out of range* (p > 1, p ≤ 0) should set
  `__se_from_clipped_pval`; *floored at the numeric bottom* should set a distinct flag with its own
  policy (default `keep`, since the SE derived from 1e-300 is still a reasonable bound).
- **Minimal** — set `pvalue.clip_low: 5e-324` so the flag fires only on true underflow.

---

## 8. The p-value floor guard is unreachable dead code

**Severity: High** · A capped −log10 p is silently flattened with `n_above = 0` and no QC line,
exactly when the user follows the documented recommendation.

**Where:** `src/postgwas/modules/harmonisation/p_values.py:553-563` and
`src/postgwas/core/statistics.py:23-27`

**Current (verified) —** `negative_log10_to_raw_p` can never return 0:
```python
    maximum_lp = -math.log10(minimum)
    value = pl.col(column).cast(pl.Float64, strict=False)
    return (
        pl.when(value > maximum_lp)
        .then(pl.lit(minimum))          # <-- returns `minimum`, never 0
        .otherwise(10.0 ** (-value))
        .alias(output_name)
    )
```
so in `p_values.py`:
```python
    # 10 ** -LP underflows to exactly 0 for very large LP; that is a p-value
    # clipped at the bottom just as surely as the cap above.
    underflow = pl.col("__raw_p") == 0        # <-- unreachable
    ...
            (high_mask.fill_null(False) | underflow).alias(CLIPPED_LOW_COLUMN),
```
The comment describes behaviour that the helper it calls specifically prevents. With
`pvalue.mlogp_max = null` — the value the policy's own help calls *"Recommended"* — `high_mask` is
`pl.lit(False)`, so LP = 400, 700, 1200 all become p = 1e-300 with nothing recorded.

**Change to:** compute the floor mask explicitly rather than inferring it from a zero that cannot
occur:
```python
    import math
    floored = pl.col("__lp_clipped") > pl.lit(-math.log10(policies.pvalue.clip_low))
    df = df.with_columns(
        [
            pl.when(floored)
            .then(pl.lit(policies.pvalue.clip_low))
            .otherwise(pl.col("__raw_p"))
            .alias(output_col),
            (low_mask.fill_null(False)).alias(CLIPPED_HIGH_COLUMN),
            (high_mask.fill_null(False) | floored.fill_null(False)).alias(CLIPPED_LOW_COLUMN),
        ]
    ).drop(["__lp_clipped", "__raw_p"])
```
**Also** add a cross-check in `policies._cross_check` requiring
`mlogp_max <= -log10(clip_low)` whenever `mlogp_max` is not null.

---

## 9. Reject accounting stops at step 13

**Severity: High** · A wrong genome build or a `chr1`-vs-`1` contig mismatch can discard up to 100%
of variants while the log reports `removed=0` and the run reports success.

**Where:** `service.py:1563, 1587, 1608, 1645`; `adapters/gwas2vcf/gwas.py:236, 244, 252, 259, 266,
274, 371, 383, 392`; `gwas2vcf_runner.py:59-74`

**Current:** `reconcile(rows_read, df.height, rejects, ...)` runs at `service.py:1563` — **before**
step 14. Steps 14, 15 and 16 then each declare:
```python
        ctx.set_rows(df.height, removed=0)
```
asserting that TSV export, gwas2vcf conversion and the whole bcftools chain remove nothing.

Meanwhile the adapter drops variants at nine points with only `logging.debug` plus an aggregate
counter — including `gwas.py:383`, where a study allele matches **neither** orientation against the
FASTA. And `run_gwas2vcf` never passes `--log`, so `main.py:92-96` defaults to `INFO` and those
debug lines are **never emitted at all**.

**Change to (three parts):**

1. **Make the adapter's losses visible.** Pass `--log DEBUG` (or a configured level) from
   `gwas2vcf_runner.py`, and have `Gwas.read_from_file` write a per-variant reject TSV using the
   reasons already registered in `rejects.py` (`invalid_position`, `non_standard_allele`,
   `reference_unmatched`, …) rather than only incrementing a counter.

2. **Replace the fabricated zeros.** The counts already exist —
   `adapter_summary.num_rows`, the adapter's `TotalVariants` / `VariantsNotRead` /
   `VariantsNotHarmonised`, and the `ORIGINAL → NORM → ID → EAF → CSQ → LIFTED` chain collected at
   `vcf_processing.py:313`. Use them:
   ```python
   ctx.set_rows(rows_before, removed=rows_before - measured_rows_after)
   ```

3. **Add a terminal assertion after step 15:**
   ```
   VCF records + VariantsNotRead + VariantsNotHarmonised == exported TSV rows
   ```
   and fail the chromosome above a configured loss fraction
   (new policy, e.g. `vcf.max_adapter_loss_fraction`, default `0.02`).

**Also fix the bare except that hides the cause** — `adapters/gwas2vcf/gwas.py:99-108`:
```python
        try:
            ...fasta fetch...
        except:
            assert 1 == 2
```
This maps a *structural* failure (contig absent from the FASTA — i.e. `chr1` vs `1`) onto the same
signal as a genuine allele mismatch. Catch `KeyError`/`ValueError` explicitly and raise a distinct
`ContigNotInReference`, so a naming mismatch fails loudly on the first variant instead of silently
discarding the chromosome.

---

## 10. `is_float32_lossy` turns a garbage variant into the top hit

**Severity: High**

**Where:** `src/postgwas/modules/harmonisation/adapters/gwas2vcf/vcf.py:22-35` and `:176-178`

**Current (verified):**
```python
    @staticmethod
    def is_float32_lossy(input_float):
        if (
            input_float == 0
            or input_float is None
            or input_float == np.inf
            or input_float == -np.inf
        ):
            return False
        out_float = np.float32(input_float)
        return out_float in [0, np.inf, 0, -np.inf]
```
Executed: `is_float32_lossy(1e-50) → True` (underflow) and `is_float32_lossy(1e50) → True`
(overflow). Both are treated identically, and then:
```python
                if Vcf.is_float32_lossy(result.se):
                    result.se = np.float64(np.finfo(np.float32).tiny).item()   # 1.175e-38
                    logging.warning(
                        f"Standard error field cannot fit into float32. Expect loss of precision for: {result.se}"
                    )
```
An SE that **overflowed** — an utterly uninformative estimate — is rewritten to the smallest
positive normal float, so `Z = ES/SE` becomes astronomically large. The message says "loss of
precision", not "value replaced".

**Change to:**
```python
    @staticmethod
    def underflows_float32(value):
        if value is None or value == 0 or not np.isfinite(value):
            return False
        return np.float32(value) == 0

    @staticmethod
    def overflows_float32(value):
        if value is None or not np.isfinite(value):
            return False
        return not np.isfinite(np.float32(value))
```
Then in `write_to_file`:
- `overflows_float32(result.se)` → **reject the variant** (count it; do not substitute).
- `underflows_float32(result.se)` → reject, or clamp with an explicit "value replaced" warning.
- `underflows_float32(result.b)` → **reject**. Currently ES below ~1e-45 is only warned about and
  then written into a `Type=Float` field as exactly `0.0`, after the `beta_zero` reject stage has
  already run.

---

## 11. A `PARTIAL` genome is reported as a success

**Severity: High** · A VCF silently missing whole chromosomes is consumed downstream as complete;
CI and batch scripts see exit 0.

**Where:** `src/postgwas/modules/harmonisation/cli.py:536-588` (`_run_validated_rows`)

**Current (verified):** `run_harmonisation_pipeline` returns `result["status"] ∈ {OK, PARTIAL}`
(set at `service.py:3552` via `_combined_dataset_status`). The CLI never reads it — grepping every
`result[...]` access in `cli.py` finds only `result["manifest"]` and `result[input_build]`:
```python
                results[row.dataset_id] = result
                _append_run_log(dataset_log, "INFO", "Dataset completed successfully.")
                ...
                print(screen_line(
                    "success", "Dataset %s completed" % row.dataset_id, indent=4,
                ))
                break
```

**Change to:**
```python
                status = str(result.get("status") or "OK").upper()
                if status != "OK":
                    message = (
                        "Dataset %s finished with status %s — the merged output is not a "
                        "complete genome. See %s."
                        % (row.dataset_id, status, result.get("manifest"))
                    )
                    _append_run_log(dataset_log, "ERROR", message)
                    _append_run_log(top_log, "ERROR", message)
                    print(screen_line("error", message, indent=4))
                    failures.append((row.dataset_id, message))
                    break
                results[row.dataset_id] = result
                _append_run_log(dataset_log, "INFO", "Dataset completed successfully.")
                ...
```
so `main()` returns non-zero. If a caller genuinely wants to accept a partial genome, gate it
behind an explicit `execution.accept_partial_genome` policy rather than on silence.

---

## 12. No memory bound anywhere in the fan-out

**Severity: High** · OOM kills, which then trigger the mis-attributed `BrokenProcessPool` retry and
re-run the whole genome.

**Where:** `service.py:874-921` (`_derive_parallelism`), `vcf_processing.py:448`,
`config/models/execution.py:10-12`

### 12a — `execution.memory_gb` is accepted and never read

Verified by grep: `memory_gb` and `usable_memory_fraction` appear **only** in the two CLI argument
tables (`cli.py:65`, `concordance/cli.py:31`) and the pydantic model. Nothing consumes them.

**Change to:** bound `workers` by memory as well as CPU. The largest chromosome partition size on
disk is already known at `service.py:2325`:
```python
    per_worker_gb = (largest_partition_bytes / 1e9) * FRAME_EXPANSION + (threads * sort_memory_mb / 1024)
    memory_workers = max(1, int((memory_gb * usable_memory_fraction) // per_worker_gb))
    workers = min(workers, memory_workers)
```
and emit the memory decision in the same `logger.decide` block as the CPU decision.

### 12b — polars' thread pool is unbounded, so the CPU budget is fictional

Under `spawn`, each worker re-imports polars, which sizes its Rayon pool from **all** visible CPUs.
`POLARS_MAX_THREADS` / `RAYON_NUM_THREADS` are set nowhere (grepped). The `threads` value computed
by `_derive_parallelism` is passed only to bcftools.

**Change to:** set them before polars is imported in the worker, via the executor initializer:
```python
def _init_worker(threads):
    import os
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    os.environ["RAYON_NUM_THREADS"] = str(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)

executor = ProcessPoolExecutor(
    max_workers=workers,
    mp_context=get_context("spawn"),
    initializer=_init_worker,
    initargs=(threads,),
)
```

### 12c — `bcftools sort -m` is given a per-thread budget for a single-threaded sort

**Current (verified)** `vcf_processing.py:448`:
```python
        sort_mem_mb = max(1, int(threads * sort_memory_mb))
```
`bcftools sort` is single-threaded and `-m` is a **total** budget. Peak sort memory across the fan
out is `workers × threads × 512 MB` ≈ `total_cpus × 512 MB` — 32 GB on a 64-core node, for sorts
that gain nothing from the multiplication.

**Change to:**
```python
        sort_mem_mb = max(1, int(sort_memory_mb))
```
and rename the policy from `vcf.sort_memory_mb_per_thread` to `vcf.sort_memory_mb`.

### 12d — the CPU budget is unset for every caller except the CLI

`execution.total_cpu_budget` defaults to `null`, and `cli._module_policy_block:219` is the only
place it is ever set. Verified at `service.py:890`: when the budget is null, `workers × threads` is
explicitly uncapped, so a direct API or pipeline caller gets `workers × 5` threads.

**Change to:** resolve it from `execution.threads` inside `_resolve_policies` (`service.py:510`) so
every entry point agrees, and reverse the `deep_merge` order in `_module_policy_block` so `--threads`
takes precedence over a run-config value (matching `CLI_OVERRIDE_PATHS` precedence used everywhere
else). While there, add the missing mapping:
```python
    configured = {
        "execution": {"total_cpu_budget": config.execution.threads},
        "logging": {
            "level": config.logging.file_level,
            "screen_level": config.logging.console_level,   # <-- currently missing
        },
    }
```
Without that second key, `--log-level` silently does nothing for the harmonisation run itself.

---

# P1 · Correctness and robustness

**13. Verify the resource build in preflight.** `resource_preflight.py:210-248` already parses every
`.fai` contig **length** and uses it only for a `length <= 0` check, then discards it. Compare those
lengths against a per-build contig-length table (new `resource_layout` entry) and fail on mismatch.
Do the same with the `##contig=<ID=…,length=…>` lines already fetched by `_vcf_info_definitions`.
Today a GRCh37 FASTA sitting in a `GRCh38/` tree passes preflight and the resulting mass REF
mismatch is **blamed on the study**.

**14. Record real provenance.** `service.py:2997-3017` records python/polars/pandas/platform but not
bcftools/htslib (the pipeline uses `--write-index` ≥1.18 and `--pair-logic` ≥1.12 features with no
minimum enforced), no checksum of any reference resource (`core/resource_preparation.sha256` exists
and is unused), and not the seed. `manifest["defaults_file"]` holds `str(dict)` of the engine config
rather than a path. Add `bcftools --version` to `require_binaries` with a minimum, `sha256` (or
size+mtime for multi-GB files) of every preflighted resource, and write the manifest atomically
(`.tmp` + `os.replace`) — every other terminal-state writer in this codebase already does.

**15. Fix the liftover gate and its tag list.** `vcf_processing.py:529-555` computes
`REJECTED / CSQ`, but variants discarded by the `exclude` SWAP filter never enter `reject_vcf` — the
YAML help states this outright. Compute on total attrition `(CSQ − LIFTED) / CSQ`. And under
`update_tags`, `--af-tags AF --es-tags ES` omits `FORMAT/EZ` and the four external `INFO/{POP}`
columns written by STEP2 — either move that annotation after liftover or extend the tag list and
sign-flip `EZ`.

**16. Resolved — reverse-complement is gated to SNVs.** Indels now test direct/swapped exact
alleles only; their representation is normalized later in the VCF workflow.

**17. Resolved — reference-AF duplicate handling is shared.** Exact allele-key/value copies keep
one by default; different AF values for the same key discard the complete key or fail by policy.
File order can no longer choose a scientific value.

**18. Add a dataset-level AF-discordance gate.** `strand.af_discordant` is written into `qc_info`
and never compared against anything (grepped). Add `strand.max_af_discordant_fraction` (fail above
~0.3). This, plus P0-A's `af_discordance_action: reject`, is what would catch a transposed
effect/other allele column — the failure mode behind Critical #2.

**19. Move numeric coercion before the mandatory-field mask.** `summary_statistics_io.py:964-977`
filters on the **raw string** frame, before `normalise_summary_statistics_values` casts at `:317`.
Cells that fail the cast (locale decimal commas, `1,234`, `NULL`) become nulls with a log line only
— no reject rows, no counter, no policy. `validate_content`'s `not_numeric` check cannot catch them
because it runs afterwards on already-`Float64` data. Coerce first, route coerced-away cells through
`RejectCollector` with a new reason `unparseable_numeric_value`, and add
`validation.on_unparseable_numeric: fail|reject|warn`.

**20. Add INFO scale detection.** `imputation_quality.py:99-104` clips everything >1.05 to exactly
1.0 with a warning, and there is no fraction gate (100% out-of-range is treated like 0.1%). A
percent-scale Rsq column is converted to uniform perfect quality and every variant then passes
`filter.info_cutoff`. Detect the 0–100 scale (median > 1.5) explicitly and add
`info.invalid_fraction_cutoff`.

**21. Recover SE row-wise, not column-wise.** `standard_error.py:153` tests `has_se` at the
**column** level, but step 06 creates the `SE` column and leaves it null wherever Z was 0. Step 08
then returns without deriving anything. A variant with a usable BETA and a usable p is rejected as
`se_null` even though `SE = |β/Φ⁻¹(1−p/2)|` is computable. Change the short-circuit to a row-level
condition: derive wherever `SE.is_null() & p.is_not_null()`.

**22. Handle `p == 1.0` / `z == 0`.** `standard_error.py:560-576`: `norm.isf(1.0/2) = 0.0` exactly,
so `se = |beta/0|` → `inf`/`NaN` → nulled → rejected as `se_null`. `validation.beta_zero = keep` is
documented as protecting β=0, p=1 but is bypassed because the variant dies on SE, not on beta. And
`pvalue.clip_high = 1.0` routes every p in (1, 1.05] into the same trap. Detect `z == 0` explicitly
and route through `validation.beta_zero`.

**23. Add a one-sided-p detector.** `se_tail` is a static policy with no data-driven check anywhere
(grepped: the string appears only in `standard_error.py` and the YAML), and the one identity that
could catch it uses `validation.z_pval_tolerance_log10 = 1.0` while a one-sided column disagrees by
exactly `log10(2) = 0.301` for every variant — so the shipped tolerance makes it undetectable by
construction. Add a study-level detector analogous to `detect_or_standard_error_scale`, and tighten
the tolerance below 0.301.

**24. Mark the Z↔p concordance check as non-independent.** `effect_validation.py:646-696` compares
`−log10(2Φ(−|Z|))` against the p it was built from whenever step 08 derived the SE — algebraically
`Z = sign(β)·z_p`, so it cannot fail. It reports "every comparable variant agrees", a **false
quality signal** on exactly the SE-less studies whose statistics are least trustworthy. The
`__se_from_clipped_pval` machinery already tracks which rows those are; report the check as
`not_independent` for them.

**25. Report `not_comparable` in the population check.** `population_frequency.py:221-222` sets
`status="not_comparable"` and appends **no warning at all** when the configured column is not in
`population_fields` — the ancestry mismatch check silently disables itself. Emit a warning.

**26. Accept a per-dataset ancestry.** `comparison_af.column: EUR` is applied to every dataset for
MAF/EAF confirmation and post-orientation AF concordance, including the conservative palindromic
AF QC. It no longer selects strand. Nothing reads a study ancestry (grepped
`ancestry|superpop|population` across `sample_sheet*` and `study_properties` — zero hits). Add an
optional `population` column to the sample sheet and select the panel column from it.

**27. Fix `chromosome_expression`.** `core/dataframes.py:39-55`: `strip_leading_zero` is
`str.replace(r"^0", "")` — one zero only, so `"007" → "07"` and `"0" → ""` rather than null — and
the function never upper-cases, while `coordinates.py:216` appends `.to_uppercase()` and
`strand.py:152` / `genome_build.py:113` do not. Use `^0+` and fold case inside the shared helper so
all four call sites agree. Also replace the fourth inline copy at `input_validation.py:1183`, which
strips the `chr` prefix unconditionally and ignores `chromosome.strip_chr_prefix`.

**28. Add an `EA == OA` reject reason.** Grepped `rejects._REASONS`, `coordinates.py`,
`input_validation.py` and `effect_validation.py`: there is no check that the two alleles differ.
With strand enabled such rows fall out as `reference_unmatched`; with it disabled they reach the
adapter and are dropped with a debug message, breaking the reconciliation in P0-B #9. Add
`monomorphic_allele` alongside `null_allele` in `coordinates.py`.

---

# P2 · Performance, structure, cleanup

**29. Stream the chromosome partitioning — the single worst memory bottleneck.**
`chromosome_partition.py:81-109` holds four whole-genome copies at once: the joined snapshot, the
snapshot `partition_by(as_dict=True)` dict, the working-frame dict, and the working frame itself.
For 20M variants that is ~20–40 GB **before a single worker starts**.
```python
    # replace both partition_by(..., as_dict=True) calls with:
    for chrom in observed_chromosomes:
        subset = df.filter(pl.col(chr_col) == chrom)
        subset.write_csv(out_path)
        del subset
```
Write the partitions compressed and read them back with the matching codec.

**30. Stop reading the input three times.** `summary_statistics_io.py:1205-1254` walks every line of
the (possibly gzipped) file in pure Python for a truncation heuristic and a cosmetic line count —
and runs even when `input.check_truncation: false`. Take the count from `df.height`, bound the
truncation check to a tail `seek`, and cache the `DelimiterDetectionResult` (delimiter resolution
currently runs up to **four** times over the same first 500 lines, each a fresh decompression).
Convert the dataset stage to `scan_csv` + one `collect()`; the concordance sub-package already
demonstrates the right pattern at `concordance/service.py:189-236`.

**31. Sample the frame for build inference.** `service.py:2040` passes the **full** study frame to
`infer_genome_build` despite the function's own docstring saying "Subset", and it then loads each
build reference fully, joins, adds two RC string columns over the joined frame and materialises four
`.unique()` frames per build. Build inference converges on a small random sample — use
`df.sample(n=min(height, 200_000), seed=…)` and count matches with `is_in`/semi-join.

**32. Filter reference panels to the working chromosome.** `allele_frequency.py:949`,
`imputation_quality.py:232` and `strand.py:125` each read a whole reference eagerly
(`pl.read_csv`, not `scan_csv`) with no chromosome filter — while `resource_paths.py:88` explicitly
invites *"one existing file that contains all chromosomes"*. With 22 workers that is 22 × the whole
panel resident and O(N_panel) work repeated 22 times. Use `scan_csv().filter(chr == chromosome)`
before `collect()`.

**33. Chain the bcftools steps.** `vcf_processing.py:336-417`: STEP0 and STEP1 correctly pipe
internally with `-Ou`, but STEP1→STEP2→STEP3 materialise `_ID`, `_EAF`, `_CSQ` as separate
bgzip+tabix files, three of which are deleted unconditionally at `:940-947`. That is three full
compressions and three index builds per chromosome on data that is thrown away. Chain them into one
`bash -c` `-Ou` pipeline, as STEP0 and STEP4 already do.

**34. Do not read a multi-GB transcript into memory in the error handler.**
`core/processes.py:186-192` reads the **entire** stdout transcript just to take the last 20 lines —
and `vcf_processing._run_bcftools_step` opens it in **append** mode, so it accumulates every stage's
stderr for the whole chromosome. A `bcftools csq` run emitting millions of GFF warnings turns a
recoverable error into an OOM inside the handler, in a worker process. Seek to the last ~64 KB.

**35. Attribute `BrokenProcessPool` correctly.** `service.py:1771-1814` labels **every** outstanding
future "worker process died" when one worker is OOM-killed, prints `_log_tail()` for chromosomes
that never started, and reruns all of them — turning one chromosome's OOM into up to 3× the entire
fan-out cost. Track which futures actually started (a start sentinel, or `future.running()` at the
break) and mark the rest *not attempted*.

**36. Wire resume.** `core/completion.py` is a complete checksummed resume framework with **zero
call sites**, while `cli.py:534` retries the entire dataset — re-reading and re-partitioning a 20M
variant file — on any `PipelineError`. Write a per-chromosome completion manifest with input/output
fingerprints and the configuration digest; skip validated chromosomes in round 1. Drop the
dataset-level retry for `PipelineError` (not a transient class; `ConfigError` is already excluded).
Either implement or remove `--resume`/`--overwrite` in `pipeline/cli.py:236-255`, which are
advertised and read by nothing.

**37. Move cleanup after the final checks.** `finalise_harmonisation_outputs` (post-merge step 3)
deletes every `chromosome_table` and `chromosome_source_snapshot`, but the raw-VCF QC assessment is
step 4 and `_resolved_harmonisation_outputs` runs after that. If step 4 fails, the inputs needed to
resume are already gone. Move the cleanup after both.

**38. Stop the worker writing to stdout.** `shared/runtime.py:143-154`
(`emit_high_visibility_warning`) writes unconditionally to `sys.stdout` and is reached from
`effect_type.py:677` **inside a chromosome worker** — breaking rule 2 of the module's own design,
which is the entire justification for the screen-text buffering. The other four call sites are
parent-side and fine. Write to stdout only when `logger is None`.

**39. Delete three dead-and-dangerous files.**
- `adapters/gwas2vcf/liftover_vcf.py` — unreferenced; hardcoded `/usr/local/bin/bcftools`; five
  `os.system(f"…")` calls interpolating user paths (shell injection, breaks on spaces); **zero exit
  checks**; a non-atomic `mv` over the output; an argparse description copied from a Mosdepth tool.
  If CrossMap fails, `mv` overwrites the real output with garbage and everything "succeeds".
- `modules/scripts/harmonisation_diagnosis.py` — 579 lines; line 14 hardcodes
  `/Users/JJOHN41/Downloads/...`; line 15 calls `pd.read_csv` while **pandas is never imported**, so
  importing it raises `NameError`; line 580 executes at import time. Superseded by `concordance/`.
- `core/execution/runtime.py` — 553 lines, every symbol unimported; the only `shell=True` in the
  tree; `print()`s from library code; five `except Exception: pass`.

**40. Resolve the remaining dead surfaces.** `core/completion.py`, `core/contracts.py` and
`core/variant_identifiers.py` have zero call sites (the first two are exactly what #36 and #41 need;
the third introduces a **fourth** variant-ID convention). `allele_frequency.resolve_palindromic_
frequency` (135 lines) is unreachable in production — its only caller passes `study_eaf_col=None` in
a branch every preceding branch has already claimed — while the real arbitration lives in
`strand.py`. The `zmaf` column is computed, maintained and never read. The per-dataset `policies`
JSON cell and `required_inputs` branch in `_resolve_policies` are emitted by nothing, and
`HarmonisationSampleSheetRow` is `extra="forbid"` so a `policies` column is rejected outright —
i.e. **the v2 sample sheet has no per-dataset policy override at all**.

**41. Declare a module contract.** `ModuleSpec` has no `requires`/`produces` fields, and `Artifact`
/ `ModuleResult` / `RunContext.publish` have zero call sites — the real handoff is untyped dicts
plus in-place mutation of `args.vcf`. Nothing validates that a `--vcf` is a harmonised GWAS-VCF, in
the expected build, from an `OK` (not `PARTIAL`) run — even though the manifest records exactly
that. Return a `ModuleResult` with `Artifact(kind="gwas_vcf", metadata={genome_build, status,
manifest})`. Meanwhile set `pipeline_target=False` for `harmonisation` so `--modules harmonisation`
is not offered by argparse and then always rejected by the planner.

**42. Consolidate duplicated logic.** One missing-value vocabulary applied case-insensitively
(currently three, and the dataset-level one is case-sensitive and lacks `NaN` — which is what
pandas/numpy/R write). One chromosome-normalisation helper (currently four). One delimiter
vocabulary (currently four, disagreeing on `pipe`; and `space`/`whitespace` both map to a single
`" "`, so PLINK `--assoc` whitespace-aligned files cannot be read by the study reader even though
the reference reader already handles them via `\s+`). One logging path (currently `PipelineLogger`,
`_announce`→stdout, ~15 bare `print()` sites, and `rich.Console`). One `Neff` implementation —
`sample_size.py:595` correctly raises on a lone case count while `adapters/gwas2vcf/gwas.py:314`
falls back to `neff = ncase`, understating balanced-design Neff by 2×, **and the adapter's version
is what reaches the VCF**. Derive `rejects.REASON_STEPS` from each module's `STEP_LABEL` — the
hand-maintained table is wrong for six reasons, so the audit table sorts removals out of pipeline
order.

**43. Emit one variant-ID convention.** `variant_identifiers.py` builds `CHROM_POS_EA_OA` (= ALT_REF
after orientation), study IDs are colon-substituted into REF_ALT order, and bcftools fills the rest
with `%CHROM_%POS_%REF_%ALT` (REF_ALT). Two of the three are allele-order-inverted relative to each
other, and anything that parses an ID back into alleles mis-keys ~half the variants. Emit
`CHROM_POS_OA_EA` throughout, and rebuild rather than retain study IDs whose parsed alleles disagree
with the harmonised EA/OA. The step also runs at position 12, **after** orientation at step 4, so a
retained ID can encode pre-orientation alleles.

**44. Small but externally visible.**
- `adapters/gwas2vcf/vcf.py:145` — remove the space in
  `##contig=<ID=…,length=…, assembly=…>`. VCF 4.x structured headers take comma-separated
  `key=value` with no surrounding whitespace; a strict parser reads the key as `" assembly"`.
- `adapters/gwas2vcf/vcf.py:137` — `if file_metadata is not None:` guards a block whose body reads
  `sample_metadata`. The guard tests the wrong variable.
- Convert the adapter's `assert`s to explicit `raise`s. `executables["python"]` is user-configurable,
  so under `python -O` **every REF-allele check becomes a no-op**.
- `adapters/gwas2vcf/param.py` uses marshmallow **2.x** APIs (`strict=True`, `.data`,
  `pass_original` 3-arg signature), removed in 3.x. Port, or pin `marshmallow<3` explicitly.
- `main.py:193-242` dereferences `fasta` after its `with` block closed it — this works only because
  pysam caches `_references`/`_lengths` at open. Move the call inside the block.
- `service.py:3388-3402` reads the gwas2vcf summary with a hardcoded `sep="\t"` against a
  configurable `gwas2vcf_input.delimiter`, and then **rewrites the file in place** as a side effect
  of computing a count. Use the configured separator; do the cleaning in
  `finalise_harmonisation_outputs` with temp+replace.
- `vcf_processing.py:95` tests plugin availability by whitespace-tokenising `bcftools plugin -l`, so
  any occurrence of the token `liftover` anywhere satisfies it. Probe
  `bcftools +liftover --version` and check the status. And `core/paths.py:88` accepts an explicit
  executable path on `is_file()` + `st_size > 0` without checking `os.X_OK`.

---

# P3 · Tests to add

These five would have caught most of P0. Each is small.

**T1 — `reconcile()` has zero tests.** Grepped `tests/`: `ReconciliationError` appears in no file.
The module's one hard conservation law could `return {}` unconditionally and the entire suite would
pass.
```python
reconcile(100, 95, collector_with_5_rejects)     # balances
reconcile(100, 96, collector_with_5_rejects)     # raises, message names "counted twice"
reconcile(100, 94, collector_with_5_rejects)     # raises, message names "vanished"
```

**T2 — the log-scale p-value conversion is never exercised.** `harmonise_p_values` is only ever
called with `decision="raw"` (`test_harmonisation_zero_p_se.py:32`, `test_harmonisation_config.py:238`).
Removing the `/ln 10` at `p_values.py:1015-1027` breaks nothing — while the module's own docstring
names this as the failure that makes "every p-value wrong by a factor of ln(10) in the exponent".
```python
# build P = -ln(p) for 1000 known p; call with decision="negln"
assert recovered_p == pytest.approx(original_p, rel=1e-12)
# same for neglog10
# detect_p_value_type returns "negln" / "neglog10" / raises AmbiguousPValueTypeError in the band
```

**T3 — the Zhu-formula test is a tautology.** `test_harmonisation_effect_from_z.py:49` asserts
`BETA == SE * 2.0`. Since `BETA = z/√den` and `SE = 1/√den`, that holds **identically** for any
denominator — it passes if `2p(1−p)` becomes `p(1−p)`, if `(N_eff + z²)` becomes `N_eff`, or if the
square root is dropped.
```python
# z=2, p=0.25, Neff=1000
assert result["SE"][0]   == pytest.approx(0.0515388, rel=1e-6)
assert result["BETA"][0] == pytest.approx(0.1030776, rel=1e-6)
```

**T4 — SE-from-p has no numeric assertion.** `test_harmonisation_zero_p_se.py:109,135` assert only
`is not None`. Changing `norm.isf(p/2)` to `norm.isf(p)` inflates every SE by ~13% at p=0.05 and
shifts every downstream Z, with no test failing.
```python
# beta = 0.2, p = 0.0455002639  (the exact two-sided tail at z=2)
assert se == pytest.approx(0.1, rel=1e-9)
# plus a se_tail=1 round-trip
```

**T5 — the adapter's `reverse_sign` is untested.** This is the **last** place effect direction can
invert, and it fires whenever the study's OA disagrees with the FASTA REF. The one real VCF test
(`test_harmonisation_gwas2vcf.py:195`) uses a matching reference and asserts only `record.id` and
the sample name.
```python
# EA=A, OA=G at a position whose FASTA base is A
assert record.ref == "A" and record.alts == ("G",)
assert record.samples[0]["ES"] == pytest.approx(-0.1)
assert record.samples[0]["EZ"] == pytest.approx(-2.0)
assert record.samples[0]["AF"] == pytest.approx(0.8)     # 1 - 0.2
assert record.samples[0]["SE"] == pytest.approx(0.05)    # unchanged
```

**Also worth fixing in the suite itself:**

- **Four tests can only pass on one machine.** `tests/data/harmonisation/manifest_v2.csv` hardcodes
  `/Users/JJOHN41/Documents/...` and `test_harmonisation_sample_sheet.py:44-48` asserts that
  **literal absolute path**. There is no skip marker anywhere in the suite (`pytest.mark` appears
  only as `parametrize`) and no `conftest.py`. Make the fixture path relative to the test file, or
  gate it behind a marker.
- **The 29 MB ADHD file is used only for existence and size checks**, never for a scientific
  regression — despite `tests/data/harmonisation/README.md` claiming it is "retained for scientific
  regression tests". Either write that regression test or drop the file.
- **Four tests assert on source text.** `test_harmonisation_logging.py:39`,
  `test_harmonisation_resource_preflight.py:335`, `test_harmonisation_screen.py:32,61` use
  `inspect.getsource()` and compare `str.index()` positions. `resource_preflight.py:338` asserts
  preflight runs before partitioning — which passes if either call sits in a dead branch or a
  comment, and breaks on a harmless rename. The critical scientific ordering (OR normalisation
  before allele swapping) is protected this way rather than behaviourally; the step-registry
  refactor in P2 #—see below—would let it be tested for real.
- **`test_harmonisation_manifest.py:127-170` patches 14 collaborators** around
  `run_harmonisation_pipeline`, so the only test touching the top-level entry point exercises zero
  science. `test_harmonisation_vcf_merge.py` writes `b"v"*200` as "VCF" files.
- **No property-based tests.** The module is full of algebraic invariants — swap∘swap = identity,
  revcomp∘revcomp = identity, `AF + swap(AF) = 1`, `Φ(Z) ↔ p` round-trip, `Z = β/SE` — and
  `hypothesis` appears nowhere.
- **A refactor that would pay for itself:** express the sixteen chromosome steps as an ordered list
  of `(number, title, operation, callable)` records and make `process_one_chromosome` a driver over
  it. Steps 14–16 collapse into one terminal "materialise" step that can be stubbed, which makes
  steps 1–13 unit-testable on an in-memory frame — and makes the ordering testable behaviourally
  instead of by string comparison.

---

# Do NOT change these

Several of these are subtle enough that a well-meaning refactor would silently break them.

1. **`strand.py` flip-vs-swap separation.** Output alleles always come from the reference
   (`:452-453`); β, Z and EAF are transformed **only** on swap (`:455-469`); SE is correctly left
   alone. A `reverse_complement` (non-swapped) match rewrites the letters without touching any
   statistic. **This is the single most important invariant in the module.**
2. **`pvalue_handler.py:34-41` uses `decimal.Decimal` with `Emin=-10000000` for LP.** Do not
   "simplify" it to numpy — the Decimal path is what preserves p-values far below float64's normal
   range.
3. **`gwas2vcf_export.py:100-132` audit-column requirements.** Because `strand_action` is written
   last and can never be empty, the adapter's `line.strip().split(delimiter)` cannot lose a trailing
   field and shift every column index. Do not relax the non-null/disjoint checks and **do not
   reorder `export_df`**.
4. **`rejects.py` `mask.fill_null(True)` + single `partition_by`.** This is what makes keep and
   reject exact complements, and it fixes a documented vanishing-row bug.
5. **`effect_from_z.py:429-431` uses `<=` / `>=` on the EAF bounds deliberately** — at EAF exactly
   0 or 1 the denominator is *exactly* zero and yields `inf`, not null. The comment explaining this
   is correct and worth keeping.
6. **`shared/variant_columns.py` forcing allele/chromosome/ID columns to `String` on first read.**
   This is what stops polars inferring an allele column of `T`/`F` as boolean, or `01` as integer.
7. **`shared/allele_join.py`** — the swapped reference frame re-keys `ea`↔`oa` *and* applies
   `one_minus` in the same `select`, direct orientation has deterministic precedence, and the join
   asserts row-count preservation. The deduplicator never picks among differing values.
8. **`--pair-logic exact` on both dbSNP and panel annotation.** A non-matching variant gets **no**
   annotation rather than a mis-oriented one.
9. **`vcf_processing.py:515-526`** — the zero-survivor invariant is correctly *not* configurable.
10. **The refusals.** `strand.py:220-244` hard-fails without an EAF/MAF decision, without an
    effect-type decision, and refuses to orient odds ratios whose SE scale is unresolved.
    `effect_type.py:857` raises rather than assuming a log-odds scale. `p_values.py:291-324` refuses
    to guess between the two log scales. `sample_size.py:595` refuses a lone case count.
    `pvalue.zero_missing_se = fail` refuses to fabricate an SE from a censored p = 0. These are the
    right failure modes and they are rarer in this field than they should be.
11. **The 243 KB defaults YAML is not bloat.** 59.7% of it is `help:` prose, median 980 characters
    per policy, **zero duplicated help strings**, every policy carries `origin`. Do not "compress"
    it.
12. **No hardcoded default sample size exists anywhere.** Every N comes from a column or an explicit
    configured value. This is the most common serious error in this class of tool and it is absent
    here — keep it that way.

---

*Nothing in this plan has been applied. All line references verified against the working tree on
2026-08-10.*
