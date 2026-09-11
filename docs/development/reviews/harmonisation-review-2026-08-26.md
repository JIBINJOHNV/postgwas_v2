# PostGWAS harmonisation — end-to-end review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/harmonisation/` (53 files, ~39,650 lines), its config model and shipped defaults, and the harmonisation documentation.
**Method:** read-only. No file was modified. Line numbers refer to the working tree as of 2026-08-25 23:35.

This is a fresh full pass, not a delta. Where it re-examines a finding from the 2026-08-25 review, the current status is stated explicitly.

---

## 0. Status of the six Critical findings from 2026-08-25

| # | Finding | Status now |
|---|---|---|
| 1 | EAF inversion undetected for a non-MAF AF column | **Fixed** — a new reference check, `assess_study_frequency_allele_alignment` (`strand.py:234-381`), now runs on exactly the `not eaf_is_maf` branch (`strand.py:721-731`) and raises rather than silently correcting. Both screening outcomes now get a reference-based test. *(But see C2 — its verdict is not consulted where it matters most.)* |
| 2 | `Neff := total N` for case-control | **Still present.** All four count combinations are now explicit and a *declared* binary trait is protected, but the undeclared case is unchanged. See **C1**. |
| 3 | Configured `pvalue.type` returned before the negative guard | **Fixed** — `assess_negative_pvalues` now runs first (`p_values.py:515-520`), and the same order holds in `study_properties.py:605-619`. |
| 4 | `info.out_of_range: clip` forcing values to 1.0 | **Fixed** — `clip` is no longer a member of the enum; default is `reject` (`harmonisation.yaml:4751-4778`). Only a bounded `(1.0, 1.05]` rounding rescale remains. |
| 5 | `POLARS_MAX_THREADS` never set | **Fixed** — exported in the parent before the pool starts, with a per-worker assertion that the limit took effect (`service.py:220-267`). This is now done correctly and needs no change. |
| 6 | `chromosome_partition.py` holding ~6 full-size copies | **Fixed** — one chromosome at a time with `del` after each write (`chromosome_partition.py:192, 208`); peak is now ~2.1× the dataset, and that 2× floor comes from upstream (see I11). |

Four of six fixed, one fixed-but-incomplete, one outstanding. The 2026-08-25 "Cleared" list was re-derived independently and still holds.

---

## CRITICAL

### C1 — A case-control study whose file carries only a total-N column gets `Neff := N`

**Where:** `sample_size.py:606-621`; `sample_sheet_generator.py:349-367`; `config/models/modules/harmonisation.py:179`; defaults `sample_size.trait_type: auto` (`harmonisation.yaml:3798`), `sample_size.controls_only: use_ncontrol_as_total` (`:3858`).

**What.** The only sample-size interface is `ncase`/`ncontrol`. The sheet generator routes a total-N column into the *control* slot, and the warning that would flag it fires only on the safe branch:

```python
controls = alternatives[0] if alternatives else None
if effective:                                    # only NEFF warns — plain N does not
    warnings.append("... review this mapping before harmonisation.")
```

`total_sample_size_column` aliases include `N`, `TOTAL_N`, `SAMPLE_SIZE` (`harmonisation.yaml:366-372`). The generator can only ever emit `trait_type: "auto"` (`Literal["auto"]`, models line 179). Step 05 then takes the `elif has_controls:` branch and writes `Neff := ncontrol` with no test that the trait is quantitative — `trait_type` is consulted only to forbid `ncase` under `quantitative` and to require both counts under `binary`.

**Why it matters.** Ncase 10,000 / Ncontrol 90,000: correct Neff = 36,000, the code writes 100,000 — a 2.78× overstatement. That value is the exported `FORMAT/NEF` *and* the `Neff` in the Z-only reconstruction `2·p·(1−p)·Neff` (`effect_from_z.py:1104-1111, 1160`), so every reconstructed BETA and SE is ~40 % too small while Z is unchanged. Downstream sample-size-weighted meta-analysis over-weights the study by 2.78×. The run reports success; only a `DECIDE` line records it.

**Fix (two small changes, neither touching correct data):**
1. `sample_sheet_generator.py:363` — change `if effective:` to `if alternatives:` so the total-N → control-count mapping carries the same review warning.
2. `sample_size.py:608` — treat undeclared as unsafe: raise when `trait_type == "auto"` and only a control count is configured, with a message telling the user to set `trait_type: quantitative` (if the lone count really is total N) or to supply a case count. Declared quantitative studies keep `Neff := N`; studies with both counts keep the harmonic formula.

Optionally, a hard cross-check: an odds-ratio effect column plus a controls-only sample size is contradictory and can be failed outright without any risk to correct input.

---

### C2 — Palindromic orientation is decided from the study EAF without requiring the alignment check to have been conclusive

**Where:** `strand.py:721-725` (the alignment check's gate) vs `strand.py:786-791` (the palindrome resolver's gate) and `strand.py:887` (`selected_pal = consensus_selected if consensus_known else frequency_selected`).

**What.** Both gates test the same three conditions — EAF column present, `not eaf_is_maf` — and the palindrome resolver never consults `frequency_alignment`. The alignment check raises only on a *conclusive* inverted diagnosis: it needs inverted r ≥ 0.80, declared r ≤ −0.80, inverted MAE ≤ 0.10, **and** a 0.10 MAE margin (`harmonisation.yaml:3510, 3538, 3560, 3581`). A column that is in fact `1 − EAF` but lands in the inconclusive band — ancestry-mismatched panel giving r ≈ −0.75, or fewer than 1,000 comparable variants on a short chromosome — passes through and then arbitrates palindromes.

When the study-wide strand consensus is *not* resolved (`consensus` is `mixed` below `consensus_threshold` 0.99, or `unresolved` below 1,000 informative variants — routine for meta-analysed sumstats mixing strands), `selected_pal` falls through to `frequency_selected`, and that suspect frequency is the sole arbiter.

**Why it matters.** Worked through the shipped thresholds: reference `REF=T ALT=A`, panel AF 0.20; study `EA=A OA=T` with a column actually holding freq(T) = 0.80. Candidates are `forward` (candidate EAF 0.80, error 0.60) and `reverse_complement_swapped` (candidate EAF 0.20, error 0.00). Error 0.00 ≤ `palindromic_af_max_difference` 0.10 and margin 0.60 ≥ `palindromic_af_min_error_margin` 0.20, so the swapped candidate wins and `strand.py:1035-1042` negates BETA and Z. Every resolvable A/T and C/G SNP — roughly 7–8 % of a typical SNV set — gets a sign-inverted effect. The post-flip AF then agrees with the panel, so the `palindromic_af_discordance` reject cannot catch it either.

Note the mechanics of the flip are correct (verified: exactly one net sign flip in all four actions, BETA/Z/EAF flipped off a single predicate, SE untouched). The defect is that the module's own safety check exists and its verdict is discarded at the one place it is decisive.

**Fix.** Add the already-computed evidence to the gate at `strand.py:786`: require `(frequency_alignment.get("declared_eaf_correlation") or 0.0) >= 0.0` — equivalently, declared MAE ≤ inverted MAE — before the study frequency may resolve palindromes. When it fails, fall through to the existing `else` branch at `:808` so those palindromes are rejected as `palindromic_orientation_unavailable` rather than oriented from a suspect column. For a correctly-oriented EAF the declared correlation is strongly positive by construction, so correct data is unaffected.

**Related, same root:** when the final frequency comes from the reference panel rather than the study (`eaf_provenance == "reference_imputed"`), the palindromic AF-discordance guard at `allele_frequency.py:1526` compares two panel values keyed on identical alleles, so it is structurally unable to fire — yet QC reports `strand_palindromic_af_comparable = N, discordant = 0`, which reads as "checked and passed". Set `comparable = 0` on that path and report `strand_qc["palindromic_resolution_counts"]` instead, so the reader can see how many palindromes rest on consensus alone.

---

## IMPORTANT

### I1 — Dataset-stage dedup deletes *both* rows of every swapped-orientation duplicate group, though chromosome step 04 exists to reconcile them

**Where:** `summary_statistics_io.py:1462-1477` (dataset pass) vs `allele_frequency.py:1466-1493` (the post-orientation pass).

The dataset pass keys on the unordered allele pair, flags groups where the *ordered* pair varies, and — under the shipped `duplicates.conflicting_action: remove_all` — rejects the whole group as `swapped_orientation_duplicate`. But the pipeline already has the machinery: after step 04 sets `EA:=ALT`, `OA:=REF`, negates BETA/Z and flips EAF, the two rows have identical ordered alleles and go through the ordinary consistency check (`post_orientation_relative_tolerance: 1e-12`) with a deterministic survivor. That later pass already handles the *harder* reverse-complement case (`A/G` vs `T/C`), which the dataset key cannot even group.

**Why:** files assembled from multiple cohorts or re-exported through several tools routinely carry the same variant in both allele orders. At shipped defaults every such variant is lost, genome-wide significant ones included, attributed to `swapped_orientation_duplicate` rather than to any scientific disagreement. The YAML rationale is right that the early pass must not *guess* — but the conclusion should be defer, not delete.

**Fix:** apply the `swapped_orientation` rejection only when `reference_aligned` is true, and exclude those groups from the early `~group_head` removal so both rows reach step 04. Irreconcilable pairs are still removed there, as `conflicting_duplicate` — the honest reason.

### I2 — An SE that overflows Float32 is rewritten to 1.175e-38 at export

**Where:** `adapters/gwas2vcf/vcf.py:180-184`; no upper-bound SE rule exists (`effect_validation.py:288-307` covers null, non-finite, and `se <= 1e-12` only).

```python
if Vcf.is_float32_lossy(result.se):
    result.se = np.float64(np.finfo(np.float32).tiny).item()
    logging.warning(f"... Expect loss of precision for: {result.se}")   # prints the fabricated value
```

`is_float32_lossy` returns True for overflow as well as underflow, so a finite SE above 3.4e38 becomes 1.175e-38 — an error of ~77 orders of magnitude, in the dangerous direction: "no information" becomes "most significant hit in the file". BETA is only warned about, never mutated; only SE is rewritten. Reachability is low (it needs a corrupt or wrongly scaled SE cell — the module's own derivations are bounded, and the `se_division_floor` guard makes underflow unreachable), but nothing anywhere stops it and there is no reject row or QC counter.

**Fix:** do not edit the vendored adapter. Add one mask beside the existing non-finite rule in `_validate_basic_effect_statistics`, under the same `validation.se_invalid` policy: reject `se.is_finite() & (se > 3.4028234663852886e38)`, with a message naming the FORMAT/SE field width. No correctly-derived SE is anywhere near that bound.

### I3 — A partly populated Z column is never completed from BETA/SE, and those rows are then rejected

**Where:** `z_score.py:76-85` and `:124-138`. The step branches on column *presence*: `if not has_existing_z:` guards the `beta/se` fill, so null cells in a supplied Z column stay null even where BETA and SE are both finite. Step 10 then removes them under `validation.z_invalid` (default `reject`).

This is the exact asymmetry the codebase already fixed for SE — `pvalue.derive_partial_missing_se` exists because "column presence is not row completeness" (`standard_error.py:209-213`).

**Fix:** in the `has_existing_z` branch, fill only null cells with the identical guarded `beta/se` expression already written at `:127-138`. Supplied Z values stay untouched. No sign question arises: step 04 has already negated BETA and the supplied Z together.

### I4 — NaN/Inf SE is "present" for recovery but "invalid" for validation

**Where:** `standard_error.py:213, 260`; `effect_from_z.py:351, 407, 671` — every "does this row need an SE?" test is `is_null()` only. Meanwhile `effect_validation.py:294-299` models `se.is_not_null() & ~se.is_finite()` as a distinct real condition with its own reject reason, so the codebase expects such cells. Nothing normalises them: `input.null_values` defaults to the exact, case-sensitive list `[NA, NAN, '', '.']`, which does not cover `NaN`, `nan`, `Inf`, `-`.

Result: a row whose SE is NaN is skipped by both recovery paths and then dropped as `se_non_finite`, even when BETA plus a signed Z, or BETA plus p, would have reconstructed it exactly. Writing `NaN` for a missing SE is routine R/pandas output.

**Fix:** change the three recovery predicates to `(col.is_null() | ~col.is_finite()).fill_null(True)`. This only widens which cells are eligible to be *filled*; every finite supplied SE is still preserved verbatim.

### I5 — The binary-trait caveat on Z-only BETA/SE is suppressed for the shipped default `trait_type`

**Where:** `effect_from_z.py:969, 1299-1311, 1332-1339` compare `trait_type` to the literal `"binary"`. The shipped default is `auto`, and the sheet generator emits only `auto`. Meanwhile `sample_size.py:593-604` computes the case-control `Neff` purely from the presence of both count columns, and `vcf_provenance._sample_size_provenance` gets it right (`if trait_type in ("auto", "binary") and has_cases and has_controls`, `:238`). The two modules disagree about the same study.

So the common case — both count columns mapped, `trait_type` left at the default — reconstructs BETA from the balanced-design Neff and exports an approximate standardized effect, while `effect_estimate_trait_type` records `"auto"`, `effect_estimate_binary_interpretation` is `None`, and no warning is emitted. A downstream user reading `FORMAT/ES` as a log odds ratio gets a systematically wrong scale with nothing in the provenance to contradict them.

**Fix:** mirror the `vcf_provenance` rule — resolve "effectively binary" as `trait_type == "binary"` or (`auto` and both count columns mapped) — and use that for the warning and the two QC keys. No number changes; only the recorded provenance and the warning.

### I6 — `_REASON_STEPS` is stale, so the two audit files disagree on which step rejected a variant

**Where:** `rejects.py:225-280` vs each module's `STEP_LABEL`. Verified mismatches include:

| Reason | `_REASON_STEPS` | Actual `STEP_LABEL` |
|---|---|---|
| `eaf_out_of_range`, `eaf_degenerate`, `af_discordant`, `palindromic_*` | `03` | `04 eaf_harmonisation` |
| `reference_unmatched`, `reference_ambiguous` | `03` | `04 strand_orientation` |
| `neff_invalid`, `sample_size_invalid` | `04` | `05 sample_size` |
| `effect_non_positive_or` | `05` | `03 effect_type` |
| `z_invalid` | `06` | `10 effect_statistics_validation` |

14 registered reasons carry a step number in `qc_summary/<dataset>_reject_reasons.tsv` that contradicts the `reject_step` column of the rejected-variant file — the two artefacts the documentation tells users to cross-reference. Separately, `eaf_unmatched_external` is registered but emitted nowhere, so its permanent zero row reads as "this check ran and removed nothing."

**Fix:** derive the map from each module's `STEP_LABEL` (single source of truth), or drop the column. Remove the unemitted reason or emit it.

### I7 — The whole-dataset retry loop catches `MemoryError` and every deterministic abort

**Where:** `cli.py:1320-1321` and `cli.py:1490-1512`. Only `ConfigError` and `ConcordanceValidationError` are special-cased, so `except Exception` catches the ambiguous-build abort, the chromosome-X Z-only block, the population-frequency inversion abort, every reconciliation failure — and `MemoryError`, which the worker layer deliberately raises specifically so the chromosome is *not* repeated under unchanged memory pressure (`service.py:1801-1806, 2405-2415`).

On a 10-million-row dataset one attempt is hours of CPU and I/O, and none of these failures is transient. The memory case directly defeats a documented safety property and risks taking the machine's other datasets with it.

**Fix:** add `except MemoryError` and `except PipelineError` handlers above the generic one, both recording FAILED and `break`ing. Leave the generic retry for genuinely transient faults (`BrokenProcessPool`, `OSError`, subprocess failures).

### I8 — Reference matching runs on un-normalised indels; `bcftools norm` runs twelve steps later

**Where:** `strand.py:643-647` (the join) and `:963` (the reject) vs `vcf_processing.py:494-499` (the only left-alignment in the pipeline, at step 16).

An indel written as `1 1000 CTT C` where the panel holds the left-aligned `1 999 ACTT AC` misses the `(chr,pos)` join entirely and is rejected as `reference_unmatched` at step 04. It never reaches the normaliser that would have made it match. The loss is systematic and allele-class-specific, and the reject reason makes it indistinguishable from a variant genuinely absent from the panel.

**Fix (minimal, honest):** do not add a FASTA read to the Polars path. Split the `reference_unmatched` count by variant class using the `snv` expression already computed at `strand.py:661-666`, expose `qc["reference_unmatched_snv"]` / `qc["reference_unmatched_indel"]`, and surface the indel figure in the step-04 QC line. If the indel share is material for a panel, normalise the *reference tables* once at staging (`external_reference_staging.py`) rather than per study.

### I9 — `reference_unmatched` conflates "coordinate absent" with "coordinate present, no compatible allele orientation"

**Where:** `strand.py:643, 933, 963`. The per-orientation candidate count is the sole discriminator, and `qc` carries no coordinate-hit counter — unlike `genome_build.py:305`, which does record `"coordinate_hits"`.

The two classes have opposite remedies (build/naming/liftover problem vs. indel representation or panel ALT coverage) and each can consume most of a study. With `unmatched_action: reject`, an operator seeing 1.2 M `reference_unmatched` cannot tell which, and the natural reaction to an allele problem is to re-check the genome build.

**Fix:** capture `coordinate_hits = joined.select(row_col).unique().height` next to the join and add it to `qc`, so `coordinate_hits − resolved − reference_ambiguous` reads as the allele-incompatible count. No retained variant changes.

### I10 — The constructed variant ID uses opposite allele order in the two VCFs one run delivers

**Where:** `variant_identifiers.py:62-74` builds `CHROM_POS_EA_OA` — and after step 04, EA is the panel ALT, so this is `CHROM_POS_ALT_REF`. Step 16 strips every ID (`bcftools annotate -x ID`, `vcf_processing.py:496-499`) and rebuilds the fallback as `+%CHROM\_%POS\_%REF\_%ALT` (`:522-524`).

Both files are promised primary outputs (`service.py:1078-1104`), so an ID-based join between them silently fails for exactly the variants dbSNP could not name — the novel and rare ones. Separately, `-x ID` is unconditional, so a study rsID absent from the configured dbSNP release is discarded rather than retained (it survives only in `FORMAT/ID`).

**Fix:** swap the two arguments in `_identifier_expr` so both delivered VCFs use the VCF-conventional `REF_ALT` order. This alters only the synthetic ID of variants that had none, only in the raw adapter VCF; no statistic, allele or coordinate changes. Retaining study rsIDs through `--columns +ID` is a separate decision.

### I11 — Performance: three changes worth most of the available win

Threading is now correct and needs no change (see §0 item 5). The remaining costs, in order:

**(a) The all-Utf8 source snapshot is pinned across five dataset steps.** `input.schema_inference_rows: 0` means every column is read as text (deliberate and scientifically right — do not change it), and `summary_statistics_io.py:1598-1599` clones that frame as `source_snapshot`. The first reject filter forks the buffers, and from then two whole-genome frames are resident until `service.py:3063`. The text frame is ~2–2.5× the typed one — roughly 6 GB against 2.7 GB at 20 M × 15 — giving the single-threaded dataset stage a ~2.2× floor, on a stage with no memory guard at all (`execution.memory_gb_per_chromosome` only sizes worker count).
*Fix:* spill it. Write the snapshot once to zstd Parquet and pass the path onward. `RejectCollector` already prefers a path (`rejects.py:337-348`), which is what per-chromosome collectors use. Only `chromosome_partition.py:166-173` needs adjusting, to `pl.scan_parquet(path).join(identifiers.lazy(), how="semi")`.

**(b) Genome-build inference holds both references and scans the fanned-out join six times.** `genome_build.py:208-218` loads both build references in one dict comprehension (both resident for the whole of `match_evidence`); the study is joined on coordinates only, so multi-allelic sites fan out 1.3–2.5×; then `:284-301` runs six separate `.unique()` passes plus three `.filter()` copies over that larger-than-study frame — ~9 extra full passes per build, ~18 for a 2-build `auto` run, each carrying two freshly allocated reverse-complement Utf8 columns.
*Fix:* loop one build at a time (the function returns only counts plus a one-column frame), and collapse `:284-301` into one `select` of `pl.col(study_row).filter(cond).n_unique()` expressions. Do **not** sample the study — those counts are user-facing QC.

**(c) The whole-genome VCF projection is written twice.** `core/vcf.py:266-284` streams `bcftools query` to a temp file, then re-reads and re-writes it in *text* mode purely to prepend a header — a full extra write, read and UTF-8 transcode of a genome-scale table, at two whole-genome call sites (`population_frequency.py:430`, `concordance/service.py:817`).
*Fix:* `run_checked_command` already supports `append_stdout` (`core/processes.py:157-160`). Write the header in `"w"` mode, then run the command with `stdout_path=destination, append_stdout=True` and delete the copy loop. Do not pass `stdout_header` — it also emits an exit-code footer that would corrupt the TSV.

Four more, each a one-to-few-line change with no result change: drop the `.to_list()` at `summary_statistics_io.py:1507` (Polars `is_in` takes a Series — currently ~200-400 k Python int allocations per pass, once per dataset *and* once per chromosome); move the `n_unique().over()` window at `:1236-1240` past the `filter` at `:1243` (it is computed over 8–30 M rows and discarded for ~99 % of them; the values are identical on the filtered frame); swap the join above the casts at `rejects.py:552-562` (currently the entire source is cast to Utf8 before the semi-join that selects a fraction of a percent of rows); and gate `_inspect_text_file` (`summary_statistics_io.py:2106-2109`) on `input.check_truncation` — it is currently unconditional and fully decompresses the input a second time, 30–60 s per dataset, for one integer used in a warning.

Also worth doing when the area is next touched: `chromosome_partition.py:156-172` does one full-frame filter *and* one hash build over the full snapshot per chromosome (24 of each); attaching the chromosome label to the snapshot once before the loop turns 24 hash joins into 1 join + 24 filters. And a whole-genome user `eaffile`/`infofile` is re-staged per *dataset* (`service.py:3078-3095`), so 40 studies sharing a panel pay for it 40 times — key the staged path on the source's content identity under a run-level directory instead.

### I12 — Error handling: nine dropped exception chains and four uninformative aborts

- **Chain dropped.** `raise X(...)` inside `except` without `from exc` at `service.py:890, 920, 933, 955, 1582, 4648, 4706`; `effect_validation.py:723`; `policies.py:653, 937`. The worst is `service.py:4706`, where a permission-denied or read-only `mkdir` becomes a `PipelineError` whose OSError shows only as incidental context — the CLI formats `type(exc).__name__, exc` into the run log, so the user loses the syscall frame. Append `from exc` at each (`from None` is the right choice at `policies.py:653`, where the ImportError adds nothing).
- **`gwas.py:459`** is separately broken: `logging.error(msg, exception_name)` passes a second positional as a %-format arg to a message with no specifier, so the logging module prints `--- Logging error ---` instead of the message. Use `%`-style args and a bare `raise`.
- **The two INFO aborts that stop a whole dataset name neither chromosome nor column** (`imputation_quality.py:866-870, 928-932`) — although nine of the eleven raises in that same function already interpolate the chromosome. On a 40-study, 23-chromosome run the user gets a bare count.
- **`"❌ Empty dataframe"`** (`summary_statistics_io.py:241-242`) is the terminal check on the primary input read and names neither the file, the resolved delimiter, nor the comment prefix — all three in scope. The most likely cause is the third: `comment_prefix="##"` consuming every data line of a file whose rows begin with `##`.
- **The vendored adapter drops variants at DEBUG level** (eight sites in `adapters/gwas2vcf/gwas.py`, launched with no `--log` so DEBUG is discarded). The user's only signal is `gwas2vcf_runner.py:134-140`: "1,204,332 exported, 1,204,109 in the VCF" with four candidate causes and no evidence. The reasons exist in `metadata` and were computed per row. Add a per-reason tally reported at WARNING after the read loop, and name the adapter transcript (already in scope as `transcript`) in the runner's error. This also restores the module's own invariant — those variants are currently in neither the output nor the reject file.
- **The delimiter cross-check disables itself silently.** `_polars_header` (`input_validation.py:339-354`) returns `None` on *any* exception, and the consumer's `if effective is not None` guard then skips the only check that the resolved separator produces the column count Polars will see. The run prints "Header check passed" and proceeds to a full parse — precisely what the three-stage design exists to avoid. A ragged-line abort from the wrong separator is the exact condition being tested and is swallowed identically to a permissions error. Return the reason and record it as a warning-severity `Problem`.

### I13 — Consolidation: two divergent chromosome normalisers, and one duplicated invariant

- **`core/dataframes.py:39-55` `chromosome_expression` is dead but discoverable, and disagrees with the real one** (`shared/variant_columns.py:45-78`) three ways: no `to_uppercase()` (so `x` and `X` stay distinct join keys — every X-chromosome variant silently lost against a reference), `r"^0"` strips one zero (`007` → `07`, not `7`), and no rename map (PLINK `23`/`26` never mapped). It is in `core/`, exported in `__all__`, and `filtering` and `formatting` already import from `core/`. **Fix:** delete both dead `core/dataframes.py` helpers and *move* the canonical `canonical_chromosome_expression` / `canonical_allele_expression` / `canonical_position_expression` there unchanged, re-exported from `shared/variant_columns.py` so the eight call sites are untouched.
- **`PALINDROMIC_SNP_PAIRS`** is a labelled biological invariant in `shared/variant_columns.py:22-25` and a hardcoded bcftools string in `filtering/sumstat_filter.py:867-872`. They agree today; a grep for the constant does not find the filtering rule. Move the constant to `core/` alongside the above and build the bcftools clause from it (byte-identical output).
- **"Count records in an indexed VCF" exists twice with different fallbacks.** `vcf_processing.py:212-239` tries `bcftools index -n` and gives up; `filtering/sumstat_filter.py:629-662` falls back to a header scan. In harmonisation, `None` is fatal (`gwas2vcf_runner.py:127-132` kills the chromosome), so on a VCF with a missing or stale index, filtering succeeds where harmonisation dies. Move the filtering implementation into `core/vcf.py` (which already owns `read_vcf_header`, `extract_vcf_table`) and keep a thin wrapper for harmonisation's warn-and-return-`None` contract.
- **The reconciliation arithmetic is written twice**, and the two disagree on severity: `rejects.py:1100-1160` raises `ReconciliationError` when the books do not balance; `sumstat_filter.py:167-176` writes `"status": "FAILED"` into a TSV and continues, producing a filtered VCF whose removal accounting is known to be wrong. `reconcile()` already accepts a plain total — call it from filtering. Note this is a deliberate behaviour change (warn → abort) and the right one: an unbalanced ledger means variants were silently lost.
- Four smaller ones: `_normalise_decision` is duplicated between `effect_type.py:114-134` and `p_values.py:419-439` (extract the *raising* form to `shared/runtime.py`; leave `sample_sheet._normalise_inferable` alone — its warn-and-downgrade is a different, deliberate contract); the concordance mismatch-reason expression is byte-identical at `concordance/analysis.py:1019-1039` and `:1624-1642`, differing only in the alias; `unused_column_name` already exists in `shared/allele_join.py` and is reimplemented inline at five sites; and two chromosome sort keys order X/Y/XY/MT differently (`service.py:808-813` vs `resource_preflight.py:49-51`) — no comparison is currently broken, but the same dataset prints two different orders.

### I14 — Both Z↔p consistency checks hardcode two-sided while `pvalue.se_tail` is configurable

`shared/statistics.py:23-32` always computes the two-sided expectation, but SE derivation honours `pvalue.se_tail` (default **2**, so shipped behaviour is coherent). Under `se_tail: 1` the pipeline is internally consistent but `validation.z_pval_concordance` compares a one-sided p against a two-sided expectation and finds every variant discordant by ~log10(2) — a flood of false warnings at the default, and a whole-chromosome removal if a site sets it to `reject`. `pvalue.se_tail` is not even in `_EFFECT_VALIDATION_POLICY_KEYS`, so the mismatch is invisible in the settings block. **Fix:** give the helper a `tail` parameter defaulting to 2 and pass `pol.get("pvalue.se_tail")` at both call sites (`effect_validation.py:735`, `effect_type.py:425-426`).

---

## OPTIONAL

- **`p` is the one statistic never derived from the others.** BETA and SE come from Z, SE from BETA and p, Z from BETA and SE — but there is no p-from-Z edge, even though the exact underflow-safe expression exists and is used for the step-10 concordance check. A study shipping BETA + SE with blank p-values in some rows loses them (`pval` is in `columns.mandatory` and in `final_check.require`). This is defensible conservatism, but it is inconsistent with `required_inputs.effect`, which explicitly documents SE as derivable. If changed: fill after step 09 into null cells only, using `LOG_P_COLUMN` provenance.
- **`pvalue.out_of_range: null` and `eaf.out_of_range: null` cannot deliver their documented "variant kept" outcome** — `final_check.require` includes both fields with `on_missing: reject`, so they are removed six steps later and attributed to the completeness gate. The authors excluded INFO from that gate deliberately and documented why; `pval` and `eaf` were not given the same treatment. Simplest fix is documentation (mirror the INFO note); the behavioural fix is to attribute the step-13 rejection to the real cause.
- **`validate_content`'s "numeric column arrived as text" diagnostic can never fire** (`input_validation.py:1183-1213`, dataset step 04) because step 03 has already cast all nine `NUMERIC_COLUMN_KEYS` to Float64/Int64. A `1,234` in a BETA column becomes a log warning, survives nine steps, and is finally rejected as `final_missing_beta` — attributing a locale problem to a completeness gate after the full fan-out. Move the check to where the coercion happens; `new_nulls` is already computed per column at `summary_statistics_io.py:361`.
- **`harmonise_chromosomes` (717 lines) inlines a 260-line reject-provenance consolidation into the retry driver** (`service.py:3566-3823`), with eight independent `raise PipelineError` sites that must each remember `results=per_chr_qc` or the caller loses the QC of every chromosome that succeeded. Extract `_consolidate_dataset_rejects(...)` with one local raise helper. No logic moves. *(The other long functions — `process_one_chromosome`, `harmonise_imputation_quality`, `derive_effect_and_standard_error_from_z` — are long because they are ordered pipelines with explicit hand-offs; splitting them would produce chains of single-use helpers. Explicitly not flagged.)*
- **`harmonise_allele_frequency` (`allele_frequency.py:1144-1780`) uses `final_df is None` as a control-flow signal** across 570 lines, and the guard at `:1464` is the honest admission that the tail's precondition cannot be read off the code. Whether `REFERENCE_AF_COLUMN` is present at `:1526` depends on which of six upstream branches ran, expressed nowhere except by column presence. Do not split the tail — return a small frozen `_EafSource(frame, column, decision_source, provenance, reference_aligned)` from the two source blocks so the precondition is explicit. Legibility fix, not a bug fix; do it only if the module is being touched anyway.
- **`sample_count_expression` (`sample_size.py:62-66`) strips `[,\s_]` unconditionally**, so a per-cohort cell `10000,20000` becomes 1,000,020,000 with no warning (the cast succeeds and `min_value` is a floor only). Contrast the EAF column, which is cast straight to Float64 so a comma-separated cell nulls with a warning. Only strip when the text matches `^\d{1,3}(,\d{3})+$`.
- **`n` and `total_n` are validated config keys that nothing reads** (`input_validation.py:797`, and the `sample_size.min_value` help text advertises them). This is exactly the misconception behind C1 — the user is told the key was accepted. Drop them, or implement them.
- **`imputation_info_column` aliases omit `R2`, `INFO_SCORE`, `IMPUTATION_QUALITY`, `IMP_QUALITY`, `MACH_R2`** (`harmonisation.yaml:378-386`). The generator writes `NA`, the sheet validates, and the run dies at chromosome step 11 after ten completed steps. Generator-only convenience list; the preflight still validates the selected column.
- **A bare `except Exception` decides whether MT variants are kept** (`coordinates.py:73-88`). Bounded today (the only pipeline call site passes `drop_mt=None`), but narrow it to `(AttributeError, KeyError)` and log which setting won.
- **Two dead branches worth a comment or a delete:** `allele_frequency.py:1457` guards on a sentinel string assigned nowhere, implying a fallthrough that does not exist (raising is the safer current behaviour and should stay); and `service.py:1911-1936` keeps a second resource-resolution path the pipeline never takes and that validates *less* than preflight — keep it as a public-API affordance but log at WARNING when it is used.
- **In the concordance audit, `np.sign(0.0)` gives Z = 0 for a zero BETA** (`concordance/analysis.py:291`), asserting a direction the data does not support. Audit-only; add `& (beta_values != 0.0)` to the usability mask.
- **`_mp_fix.py` forces `mp.set_start_method("fork", force=True)`** and is imported nowhere in the snapshot. If `src/postgwas/__init__.py` imports it, any pool created without an explicit context would fork with Polars' Rayon pool live, defeating the `POLARS_MAX_THREADS` discipline. Confirm and delete.
- **Dtypes:** positions are Int64 where Int32 is exact for human coordinates (~80 MB at 20 M rows, and 4-byte instead of 8-byte join keys); chromosome is Utf8 where `pl.Enum(chromosome.allowed)` would cut ~120 MB — prefer Enum over Categorical, since the column is compared with `==`/`is_in` across frames and Categorical comparison depends on the global string cache. **BETA, SE, Z, P, EAF and INFO are Float64 and must stay Float64** — `p_values.py` goes to considerable lengths to preserve sub-normal p-values, and narrowing any statistic would destroy that.

---

## Verified correct — do not re-flag

Re-derived independently this pass:

- **Exactly one net sign flip in all four strand actions.** `forward` → none; `forward_swapped` → one; `reverse_complement` → none (relabelling the strand does not change which physical allele is the effect allele); `reverse_complement_swapped` → exactly one. BETA, Z and EAF flip off a single `swap` predicate (`strand.py:1030-1045`); SE, INFO, N and p correctly untouched.
- **The 1−EAF flip happens exactly once**, in one of two mutually exclusive branches keyed on `eaf_is_maf` (`strand.py:1043-1045` vs `allele_frequency.py:1364-1373`); palindromic candidates are scored on a separate working column, so no double flip.
- **OR → log(OR) is ordered before orientation**, and `strand.py:581-589` raises if raw odds ratios reach it. `SE[log OR] = SE[OR]/OR` is evaluated against the pre-transform OR in a single `with_columns`.
- **METAL and Zhu et al. algebra** (`effect_from_z.py:1105-1112, 1154-1261`): the `2·p·(1−p)` factor is present, `Neff` (not raw N) is used, `zhu_2016` correctly adds `+ Z²`, and `BETA/SE` returns exactly the input Z.
- **Full log-space p ↔ Z round trip** via `ndtri_exp` / `log_ndtr`: p = 1e-320 → |Z| = 38.2872 → 320.000 exactly. The exact source token is preserved before any numeric cast, so sub-normal p-values reach the adapter intact. `p = 0` becomes `clip_low` with a marker, never an infinite Z; `p = 1` is not clamped.
- **Neff formula** `4 / (1/Ncase + 1/Ncontrol)` as written, guarded against zero, and identical in the adapter.
- **Reference-panel AF confirmation** excludes palindromes, requires 1,000 non-palindromic matches, refuses when the panel's own effect alleles are predominantly minor, and requires a 0.02 separation — it never guesses.
- **Duplicates are resolved twice and the second pass is after orientation**, with both alleles in the key (so distinct ALTs at one coordinate are never collapsed) and a retention order that `policies.py:1366-1382` *enforces* to end in `input_order`, giving a strict total order regardless of sort stability.
- **Multi-allelic sites are matched per-ALT**, not rejected wholesale; genuine ambiguity is `reference_ambiguous`, distinct from `reference_unmatched`.
- **Rejected rows are recorded before any transformation** and restored from an immutable snapshot by stable row ID, so reject files show original study values. Row preservation is asserted, not assumed (`reconcile()` hard-asserts rows-in − rejects = rows-out per chromosome).
- **Non-SNV alleles are never complemented**; the complement is applied to the strings via `replace_many` + `str.reverse` (verified `AC → GT`, not `AA`), and output alleles come from the reference.
- **Allele-independent annotations are not inverted** — external INFO joins with `swapped_value="same"`; only frequencies use `"one_minus"`.
- **Palindrome band `[0.40, 0.60]` inclusive**, cross-validated to straddle 0.5; consensus and frequency evidence cross-check each other and conflicts are rejected, never silently resolved.
- **Thread discipline** (`POLARS_MAX_THREADS` + spawn + per-worker assertion) and **`chromosome_partition` memory** are now correct.

---

# Documentation review — harmonisation only

The step *sequences* are accurate: the 8 / 16 / 5 step lists in `processing-order.md:71-101`, `README.md:251-297` and `README.md:337-343` match `service.py` exactly, including the two unnumbered gates. **No case was found where a doc presents steps in an order the code does not follow.** The defects are inside step descriptions, in what is omitted, and in a handful of wrong controls.

One premise from the previous review needs correcting: `README.md:479-760` *does* contain a complete generated policy registry — all 161 policies, 0 wrong defaults. What does not exist is a cross-cutting **disposition** view (§F below).

## Missing from the docs

**[Critical] The output is intersected with the `default_eaf` panel.** `strand.py:643` is an `inner` join on chromosome+position, `require_default_eaf` is `True` for every chromosome, and `strand.unmatched_action` defaults to `reject`. This is the single largest variant-loss mechanism in the pipeline and it is never stated as such — panel density and correct build, not strand settings, set the ceiling on output size. A user seeing a 40 % drop will investigate the settings the docs *do* explain. Add to processing-order step 4, and point at the `reference_unmatched` row of the reject report as the first thing to check.

**[Critical] `eaf.degenerate: reject` drops every variant with EAF exactly 0 or 1.** The step-4 walkthrough ends at the AF-concordance QC; the rejection that runs after it (`allele_frequency.py:1683-1710`) is absent. Monomorphic markers — common in imputed panels and in studies that round EAF — are removed silently. `configuration.md:404` mentions only `eaf.degenerate: keep`, in passing, so the shipped `reject` appears nowhere.

**[Critical] A reported `p = 0` is replaced by 1e-300, and `1 < p ≤ 1.05` is set to 1.0.** `processing-order.md:517-519` says only "applies `pvalue.out_of_range` to invalid values". Two values are *overwritten*, and 1e-300 propagates into the derived SE, Z and `FORMAT/LP`. The page never gives the default action, the substituted values, or the fact that the same `clip` setting still *rejects* non-finite, negative and >1.05 values.

**[Important] Dataset step 8 makes the INFO-vs-MaCH-Rsq decision and can abort the dataset there** — when >0.1 % of finite INFO values exceed 2.0, or the column has no finite value at all (`imputation_quality.py:186-224`). A common failure mode for a mis-mapped INFO column, and step 8's description covers only partitioning and parallelism.

**[Important] >1 % missing sample size stops the dataset**; below that, those variants are silently removed at chromosome step 5. `processing-order.md:299` says only "too much missing sample size stops the dataset" — the 1 % threshold appears in no wiki page.

**[Important] Standard-INFO values in `(1.0, 1.05]` are rewritten to 1.0** (`imputation_quality.py:810-830`). Step 11 says only "follow the resolved policy".

**[Important] The exported effect and other alleles are the panel's ALT and REF, not the study's spelling.** `processing-order.md:647` describes reverse complement as a per-row string operation; in fact every retained variant's EA/OA are replaced wholesale (`strand.py:1029-1033`). One sentence explains why exported alleles differ from the input even for `forward` rows.

**[Important] There are two independent retry mechanisms.** `execution.max_retry_rounds` retries failed *chromosomes* inside one run; `execution.retries` re-runs the *entire* dataset pipeline including the full input read. Stage 4 documents only the first, so a user will not expect duplicated dataset logs.

**[Optional]** `Neff` is rounded to a whole number (matching the integer `FORMAT/NEF`); and `--zero-p-se-action` — the only CLI escape from the shipped `pvalue.zero_missing_se: fail`, and the flag the error message at `standard_error.py:233` tells users to reach for — is documented nowhere.

## Described but not done

**[Important] `input.strip_double_hash_lines` help describes a file copy that does not happen.** The YAML says the file "is copied without those lines", that a `.gz` input "is written out UNCOMPRESSED next to the original", and that "the input folder must be writable and must have room". The code passes `##` to Polars as `comment_prefix` (`summary_statistics_io.py:133, 166-173`). No copy exists anywhere. This text is what `postgwas config export --style full` prints, and the README tells users to read it — it states a disk-space and write-permission requirement that does not exist.

**[Important] Step 12 does not use `vcf_processing.missing_id_format`.** `README.md:293` names it as the control; `variant_identifiers.py:61-71` builds a fixed string with no policy lookup, and `missing_id_format` is read only at `vcf_processing.py:522-524`. A user who changes it will see step 12 ignore it — and the two fallbacks use opposite allele order (see I10).

**[Important] `configuration.md:345-346` names the wrong skipped step.** With no effect column, chromosome step **03** (effect type) is skipped, not step 05 — step 05 is sample size and must always run, because step 06 needs `Neff`. A user matching the sentence against a log looks at the wrong `SKIP` line.

**[Optional]** `field_lifecycle` in the YAML says the constructed ID is `chr:pos:ea:oa` (the code uses underscores and explicitly converts colons *to* underscores) and that SE-from-p happens at "chromosome step 09" (it is step 08; step 09 is Z). `pvalue.se_tail`'s help repeats the step-09 error. Both strings are exported to users.

**[Optional]** `sample-sheet.md:8-16` lists `trait_type` and `delimiter` under "Required" although both default to `auto` — and line 118 of the same page says so correctly.

## Wrong order

None found at the step-sequence level. Two sub-step notes:

- `processing-order.md:295-298` reverses the two operations inside dataset step 4 — `prepare_missing_sample_sizes` runs *before* `validate_content` (`service.py:2628-2636`). It matters only because the sample-size abort fires first, so a dataset with both problems reports the sample-size failure.
- Chromosome step 4 resolves palindromic orientation before the EAF range gate. **This is not a defect** — `strand.py:425-440` independently requires a finite frequency in `[0,1]` before a candidate is scored. Recorded so it is not mistaken for one.

## Terminology, options, outputs, failures

- **`CONCORDANCE_FAILED` appears in no doc** (`cli.py:1474-1476`), yet it is exactly what a `--validate` user sees when harmonisation succeeded but the audit did not — and the VCFs exist. Without documentation it reads as total failure. The documented status set also omits `PREFLIGHT_FAILED`, `RUNNING` and `INTERRUPTED`.
- **The README output tree omits one of the three promised primary VCFs.** `_resolved_harmonisation_outputs` (`service.py:1078-1104`) requires the raw `{dataset_id}_gwas2vcf_{build}_merged.vcf.gz` as well as both merged build VCFs, and raises if any is missing. `outputs-and-qc.md:11-16` gets this right; the README tree contradicts it, so a user tidying outputs may delete a file the command validates.
- **The reject-reason report's `step` column contradicts the rejected-variant file's `reject_step`** for 14 reasons — this is the code defect I6, not a doc gap. Once fixed, say that the report carries the owning step and a plain-language description.
- **Three routinely written outputs are not listed** in `outputs-and-qc.md`: `{dataset}_field_completeness.tsv`, `{dataset}_{build}_population_frequency_qc.json`, and `{dataset}_post_orientation_duplicates_chr{N}.tsv` (written unconditionally per chromosome, even when empty).
- **`NEF` and `NCO` are missing from the GWAS-VCF field list** (`README.md:1179-1181` says "`ES`, `SE`, `LP`, `AF`, `SI`, `SS`, `NC`, `EZ` and related"). `NEF` carries the effective sample size the whole `Neff` discussion is about.
- **"clears temporary record IDs"** (`processing-order.md:595`) understates `annotate -x ID`, which strips *every* ID including genuine study rsIDs — they survive only in `FORMAT/ID`.
- **`eaf_qc/` records are post-orientation values**, not values as supplied: for a swapped row a supplied 1.5 is recorded as −0.5. The as-supplied values are in the rejected-variant file.
- **The sample-sheet preflight failure list omits duplicate input files**, which `sample_sheet.py:406-411` rejects.

## The one addition worth making

A **"What the shipped defaults do to your data"** section in `processing-order.md`, after "What can happen to one input row?". The registry in the README has every value, but spread across 20 collapsed `<details>` blocks grouped by policy family, each disposition sitting in column 2 beside ~145 thresholds and column names. Nothing lets a reader answer the two questions that matter: *which shipped defaults silently drop my variants, and which stop my run.* Grepping `docs/` for `max_missing_fraction`, `af_tolerance`, `eaf.degenerate`, `pvalue.out_of_range`, `info.on_missing` or `final_check.on_missing` returns essentially nothing.

Every value below was read from `harmonisation.yaml` and confirmed against the code site named.

### Removes variants — no configuration needed to trigger it

| Config key | Default | Consequence | Code |
|---|---|---|---|
| `strand.unmatched_action` | `reject` | **Largest single filter.** Chromosome+position absent from the `default_eaf` panel, or no REF/ALT orientation matches → `reference_unmatched`. Output ≈ study ∩ panel. | `strand.py:963` |
| `strand.ambiguous_action` | `reject` | More than one valid orientation, or a palindrome with neither strong consensus nor decisive EAF. | `strand.py:968-1027` |
| `eaf.out_of_range` | `reject` | Non-finite EAF, or EAF outside `[0,1]`. | `allele_frequency.py:349` |
| `eaf.degenerate` | `reject` | **EAF exactly 0 or 1** — monomorphic markers do not survive. | `allele_frequency.py:1698` |
| `strand.palindromic_af_discordance_action` (tol `0.2`) | `reject` | Palindrome whose oriented EAF differs from panel AF by >0.2. | `allele_frequency.py:1624` |
| `effect.or_non_positive` | `reject` | OR ≤ 0, before `ln(OR)`. | `effect_type.py:806-814` |
| `sample_size.missing_action` | `remove` | Missing case/control/total count (only reached if ≤1 % of the study). | `sample_size.py:504` |
| `effect_from_z.beta_z_sign_mismatch` | `reject` | BETA and signed Z with opposite signs, where SE must be recovered. | `effect_from_z.py:331-440` |
| `effect_from_z.low_effective_variance_action` (min `1.0`) | `reject` | Z-only reconstruction: `2·EAF(1−EAF)·Neff < 1`. | `effect_from_z.py:1099-1153` |
| `validation.se_invalid` (floor `1e-12`) | `reject` | SE missing, non-finite, or ≤ 1e-12. | `effect_validation.py:277-303` |
| `validation.beta_invalid` / `validation.z_invalid` | `reject` | BETA or Z missing or non-finite. | `effect_validation.py:311-341` |
| `validation.se_from_clipped_pval` | `reject` | SE derived from a clipped p — how `p = 0` rows usually leave. | `effect_validation.py:449-479` |
| `info.out_of_range` | `reject` | INFO below 0 or above the resolved maximum (1.05 standard / 2.0 MaCH Rsq). | `imputation_quality.py:886-903` |
| `final_check.on_missing` (require `chr,pos,eaf,beta,se,zscore,pval`) | `reject` | **Last gate** — most often EAF or SE. `info`, `n`, `snp` are not required. | `effect_validation.py:965-971` |
| `duplicates.conflicting_action` | `remove_all` | Duplicate groups with disagreeing statistics, and any group carrying both EA/OA orders, dropped **entirely** — no survivor. | `summary_statistics_io.py:1855-1878` |
| `external_reference.non_identical_duplicate_action` | `discard_all` | A reference key with disagreeing values contributes nothing; the study row behaves as unmatched. | `strand.py:159-216` |
| `chromosome.drop_mt` / `allowed_after_split` | `true` / `1–22, X` | Y, MT, scaffolds and alt contigs → `unsupported_chromosome`. | `coordinates.py:93-104, 300` |
| `vcf.liftover_swap` | `exclude` | Target-build VCF only: records the plugin marked swapped or reference-added are removed. The YAML's own `recommendation:` for this key is `keep`. | `vcf_processing.py:621-627` |

### Overwrites a value rather than dropping the variant

| Config key | Default | Consequence | Code |
|---|---|---|---|
| `pvalue.out_of_range` + `clip_low` | `clip` + `1e-300` | Reported `p = 0` → **1e-300**; propagates into derived SE, Z and `FORMAT/LP`. | `p_values.py:1024-1031` |
| `pvalue.clip_high` + `tolerance_above_one` | `1.0` + `1.05` | `1 < p ≤ 1.05` → exactly 1.0; `p > 1.05` rejected. | `p_values.py:1024-1031` |
| `info.clip_max` + `clip_tolerance` | `1.0` + `1.05` | Standard-INFO in `(1.0, 1.05]` → 1.0. | `imputation_quality.py:810-830` |
| `info.multi_value_aggregation` | `median` | A delimited INFO list collapses to its unweighted median. | `summary_statistics_io.py:372` |
| `sample_size.controls_only` | `use_ncontrol_as_total` | With no case count, the control count becomes `Neff` outright — **see C1**. | `sample_size.py:606-621` |

### Keeps the variant

| Config key | Default | Consequence | Code |
|---|---|---|---|
| `validation.beta_zero` | `keep` | A finite BETA of exactly 0 is a legitimate null result (`Z = 0`). | `effect_validation.py:319-330` |
| `info.on_missing` | `keep` | No imputation quality → still exported; INFO is not required at export. | `imputation_quality.py:955-963` |

### Warns only

| Config key | Default | Consequence | Code |
|---|---|---|---|
| `strand.af_discordance_action` (tol `0.2`) | `warn` | A **non**-palindromic EAF more than 0.2 from panel AF is reported and kept — unlike its palindromic twin above. | `allele_frequency.py:1651-1665` |
| `validation.beta_se_z_concordance` (`0.01` abs + rel) | `warn` | `Z ≠ BETA/SE` reported, not removed. | `effect_validation.py:625-632` |
| `validation.z_pval_concordance` (tol `1.0` log10) | `warn` | Z and P disagreeing by up to a **factor of ten** pass silently. | `effect_validation.py:487-497` |
| `validation.declaration_mismatch_action` | `warn` | An explicit sheet `effect_type` / `p_value_type` / `delimiter` that the detector disputes **wins**. | `study_properties.py:604-612` |
| `population_frequency_qc.on_error` | `warn` | An operational failure does not fail the dataset — but a diagnosed inversion still does. | `service.py:4258-4270` |

### Stops the dataset or chromosome

| Config key | Default | Consequence | Code |
|---|---|---|---|
| `effect_from_z.x_chromosome_z_only_action` | `fail` | X present with Z but neither BETA nor SE. | `service.py:2679-2698` |
| `pvalue.zero_missing_se` | `fail` | Literal `p = 0` with no supplied or Z-derived SE. Override with `--zero-p-se-action`. | `standard_error.py:224-238` |
| `sample_size.cases_only` | `fail` | A case count without a control count cannot give `Neff`. | `sample_size.py:624-635` |
| `sample_size.max_missing_fraction` | `0.01` | **>1 %** missing counts stops the dataset at step 4. | `sample_size.py:200-208` |
| `effect_from_z.max_beta_z_sign_mismatch_fraction` | `0.01` | **>1 %** whole-study BETA/Z sign discordance. | `effect_from_z.py:184-330` |
| `pvalue.max_negative_fraction` | `0.001` | **>0.1 %** negative usable p-values. | `p_values.py:186-205` |
| `info.maximum_invalid_fraction` (`mach_rsq_max` `2.0`) | `0.001` | **>0.1 %** of finite INFO values above 2.0, at step 8. | `imputation_quality.py:208-224` |
| `validation.max_reject_fraction` / `warn_reject_fraction` | `0.95` / `0.2` | Step 10 alone removing >95 % fails the chromosome; >20 % warns. | `effect_validation.py:522-546` |
| `validation.on_ambiguous_column_mapping` / `on_duplicate_header` | `fail` | Two config keys on one column, or a repeated header name. | `input_validation.py:1046` |
| `build.mode` | `auto` | An `Ambiguous` build verdict stops the dataset. | `service.py:2820-2826` |
| `execution.fail_dataset_on_chr_error` | `true` | One chromosome still failing after retries fails the dataset; `false` accepts `PARTIAL`. | `service.py:3915-3926` |

---

## Suggested order of work

1. **C1** and **C2** — both are small, gated changes that cannot affect correctly-handled data.
2. **I6** (`_REASON_STEPS`) and **I2** (SE upper bound) — one-table and one-mask changes.
3. **I3**, **I4**, **I5** — three small predicate/condition changes that recover variants or fix provenance.
4. **I7** — two `except` clauses in the CLI.
5. **I11(c)**, then the four one-line performance changes; then **I11(a)** and **(b)**, the two structural ones.
6. **I1** — the deduplication deferral, which is a judgement call against a documented decision; read the counter-argument in that entry first.
7. Documentation: the three Critical omissions and the disposition table, then the "described but not done" corrections (the `strip_double_hash_lines` help is the most misleading).

Nothing recommended here changes the value of a correctly harmonised variant.
