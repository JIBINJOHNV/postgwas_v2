# PostGWAS `formatter` — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/formatting/` (`service.py` 1296, `resume.py` 871, `table.py` 592, eight exporters 671, `cli.py` 196, `contracts.py` 125, `ldsc_reference.py` 204, `reference_identifiers.py` 134), `config/models/modules/formatting.py`, `config/defaults/modules/formatting.yaml`, `docs/wiki/modules/formatting.md`, `tests/test_formatting.py`, and the upstream contract with `harmonisation`.
**Method:** read-only on the tree — no file was modified. Each export's column mapping was checked against the target tool's documented input format. The performance claim in **I1** is measured on a 5,000,000-row frame using the module's own identifier expression and duplicate policy. Line numbers refer to the working tree as of 2026-08-26.

**Summary.** Every column mapping I checked is scientifically correct — GCTA `.ma` positional order, FINEMAP's `.z` field order with a true MAF transform, LDSC's A1-as-effect-allele, MiXeR's `SNP CHR BP A1 A2 N Z`, and MAGMA's two files built from one shared variant set. The module also does something no other module here does: `contracts.py` records, per target, what the frequency column means and cites the tool's own documentation. There are no Critical findings.

What the review turns up instead is **one measured inefficiency that dominates the module's runtime** (the same identifier resolution is recomputed once per target, and with the shipped defaults every target's arguments are identical), one **decision rule worth tightening** (`infer_study_design` treats a single non-null case count as proof of a case-control study), and a structural observation: 2,167 lines of orchestration and resume machinery wrap 671 lines of actual export logic.

---

## IMPORTANT

### I1 — The same variant-identifier resolution is recomputed once per target; with the shipped defaults it is identical every time

**Where:** `service.py:1046-1103` (the per-target loop), `table.py:122-201` (`select_variant_identifiers`).

Inside the target loop, every target calls `select_variant_identifiers` on the **full** frame. That call does a regex extraction over every row, a filter, an `is_duplicated()` scan, an `n_unique()` on the duplicated subset, and then `resolve_duplicate_identifiers` (another group-by):

```python
identifier_type   = module.variant_identifiers.target_types.get(target, module.variant_identifiers.default_type)
duplicate_policy  = module.variant_identifiers.target_duplicate_policies.get(target, module.variant_identifiers.default_duplicate_policy)
...
target_frame, identifier_qc = select_variant_identifiers(frame, module, identifier_type, duplicate_policy=duplicate_policy)
```

The shipped configuration sets `default_type: rsid`, `default_duplicate_policy: exclude_all`, and leaves **both** `target_types` and `target_duplicate_policies` empty (`formatting.yaml:115-123`). So for a run selecting all seven targets in `format_order`, the arguments are byte-identical seven times and so is the result.

Measured on a 5,000,000-row frame with the module's own `rsid_extraction_pattern` and the `exclude_all` policy:

| | wall time |
|---|---|
| one identifier resolution | **8.61 s** |
| seven, as the shipped defaults perform | **48.41 s** |
| avoidable | **39.80 s (82%)** |

The fix is a memo keyed on `(identifier_type, duplicate_policy)`; the LDSC-with-`--merge-alleles` path (`service.py:1063-1097`) is genuinely different and stays outside the memo. On the shipped configuration this collapses seven full-frame passes into one, and it remains correct if a user later gives one target its own identifier convention — the key changes and the memo misses.

### I2 — `infer_study_design` treats a single non-null case count as proof of a case-control study

**Where:** `table.py:518-550`.

```python
cases = int(counts["cases"])          # number of variants with a non-null NC
...
trait_type = "quantitative" if cases == 0 else "binary"
```

Two consequences, in opposite directions.

**One stray value flips the study.** If three variants out of ten million carry an `NC` value — a partially populated source column, a merge artefact — the study is declared binary. LDSC then writes `N_CAS`/`N_CON`, and `complete_rows` requires both to be non-null and positive, so the export retains **only those three variants**. To the module's credit this does not pass silently: `service.py:1017-1030` emits a warning naming exactly this condition, and the per-target QC step reports the excluded count with a reason. But the decision itself still hinges on a single row. A proportion test — or a policy of failing when case counts are present for some but not most variants — would match the care taken everywhere else in this module.

**The converse is a live cross-module risk.** The harmonisation review's still-open Critical (a case-control study whose file supplies only a total-N column has that N routed into the *control* slot at `sample_sheet_generator.py:362`, with no `NC` written) lands here as `cases == 0` → `trait_type = "quantitative"`. LDSC then receives `N` = total N with no case/control split, so no liability-scale conversion is possible and the reported heritability is on the observed scale for a binary trait, labelled quantitative. The formatter behaves correctly given its input — the point is that it is the last place in the chain that could notice, and by design it cannot: `infer_study_design`'s docstring says "never metadata", and there is no configuration key to assert the trait type. That is a defensible choice, but it means the harmonisation gap propagates all the way to a wrong h² with nothing in between to challenge it.

Worth considering an explicit optional `study_design.trait_type` override that, when set, is checked *against* the inferred value and raises on disagreement — keeping inference as the default while giving a user a way to catch exactly this.

### I3 — 2,167 lines of orchestration and resume wrap 671 lines of export logic

`service.py` (1,296) and `resume.py` (871) together are 3.2× the eight exporters combined. `resume.py` does build on `core.completion` rather than duplicating it, and much of its bulk is genuine formatter-specific work — per-target output paths, rebasing recorded artifact paths when a run directory is copied, matching historical configurations across schema changes. But two functions are large enough to be worth splitting on their own: `_rebase_manifest_results` (`:438-582`, 144 lines) and `_validate_manifest` (`:655-820`, 165 lines).

This is an observation rather than a defect — the test suite covers the resume paths thoroughly (twelve of the forty tests in `test_formatting.py` are resume cases). It is worth flagging only because the ratio makes the module look far more complex than its science is, and future scientific changes have to be made through that layer.

---

## OPTIONAL

**O1 — PRED-LD writes a file for every configured chromosome, including ones with no variants.** `export_pred_ld` loops over `config.chromosomes` and calls `write_dataframe_table` on `output.filter(CHROM == chromosome)` (`pred_ld.py:80-93`); `write_dataframe_table` writes unconditionally (`core/io/tables.py:192-207`), so an absent chromosome yields a header-only file, and that path is still returned in `files`. The imputation module consumes that list. Worth confirming PRED-LD tolerates an empty input, or filtering the loop to chromosomes actually present.

**O2 — the total-N versus effective-N split across targets is deliberate but unexplained.** MAGMA and GCTA receive `N` (total, `[%SS]`); SuSiE, FINEMAP and MiXeR receive `NEFF` (`[%NEF]`); LDSC receives case/control counts. `sample_size_reporting.target_notes` (`formatting.yaml:93-111`) documents *what* each target gets, and the screen table renders it — but not *why* MAGMA gets total N for an unbalanced case-control study, where effective N is the more common recommendation for meta-analysed summary statistics. One rationale line in the config comment would close it; if the answer is "MAGMA's manual says sample size", say that.

**O3 — `minimum_p_value: 1.0e-300` is a policy clamp, not a numerical necessity.** `10**(-LP)` remains representable to about 1e-308, so the bound is a chosen export contract. It is applied and counted honestly — `count_bounded_negative_log10_values` reports `p_values_bounded` per target — and it matters only for MAGMA and GCTA, since LDSC's `munge_sumstats` works from `Z`. Worth one sentence in the docs saying which targets the clamp can actually influence.

**O4 — `infer_study_design` scans the frame independently in `export_ldsc` and `export_mixer` when `study_design` is not supplied** (`ldsc.py:24-29`, `mixer.py:23-28`). `service.py:996` computes it once and passes it, so the fallback is dead in the pipeline path — but it is live for anyone calling an exporter directly, and two exporters computing the same study design from the same frame is the kind of thing that later drifts.

**O5 — `select_variant_identifiers` computes `duplicate_groups` via a second `n_unique()` pass** over the duplicated subset (`table.py:169-176`) purely for a QC field. It is reported and worth keeping, but it can come out of the same aggregation that already identifies the duplicates.

---

## Verified correct — do not re-flag

**Every export's column contract.** Checked against each tool's documented input format:

- **GCTA `.ma`** — `SNP A1 A2 freq BETA SE P N` with `A1 = ALT` (the effect allele) and `freq = EAF` (the effect-allele frequency). GCTA reads this file positionally and ignores the header, so the header spelling `BETA/SE/P` rather than `b/se/p` is cosmetic and the *order* is what matters — and it is right.
- **FINEMAP `.z`** — `rsid chromosome position allele1 allele2 maf beta se`, with `allele1 = ALT` (FINEMAP's effect allele) and `maf` produced from EAF by `effect_frequency_to_minor_frequency` = `min(f, 1−f)` (`table.py:508-510`). A true MAF, not the effect-allele frequency.
- **LDSC** — `A1 = ALT`, `A2 = REF`, matching LDSC's A1-is-effect convention; binary studies get `N_CAS`/`N_CON`, quantitative studies get `N` from `N_CONTROL` (which harmonisation defines as total N for quantitative traits).
- **MiXeR** — `SNP CHR BP A1 A2 N Z` with `N = NEFF`.
- **MAGMA** — the SNP-location and p-value files are both derived from the *same* `usable` set (`magma.py:20-38`), so the annotation and the association file cannot disagree about which variants exist.

**Z is always available.** A natural worry is that `Z` maps only to `[%EZ]`, and the GWAS-VCF spec defines EZ as "Z-score provided if it was used to derive the EFFECT and SE fields" — which would leave SuSiE, LDSC and MiXeR with an empty required column for the common beta/SE study. It does not: harmonisation derives `Z = BETA / SE` when the study supplies no Z (`harmonisation/z_score.py:126-147`, with an `se_division_floor` guard and nulls for unusable rows), assigns it to `imp_z_col`, and `service.py:623` maps that into the VCF writer's Z field. The formatter correctly does no derivation of its own.

**MiXeR's median-N rule is applied where MiXeR specifies it.** The median effective sample size is computed on the study-wide frame *before* the INFO and variant-type exclusions, with a comment saying exactly that (`mixer.py:51-59`), and only then is the `minimum_sample_size_fraction` cut applied. Computing the median after those exclusions would shift the threshold.

**LDSC sample prevalence.** Median (configurable to mean) of the per-variant `NC / (NC + NCO)`, with the minimum, maximum, and the case- and control-count ranges all reported alongside it (`ldsc.py:57-97`). The config comment states the rule that matters — "Variant counts must never be summed" (`formatting.yaml:144`).

**The p-value bound is counted, not hidden.** `negative_log10_to_raw_p` clamps rather than emitting 0, and `count_bounded_negative_log10_values` reports how many values hit the bound for every target that transforms LP.

**One VCF read.** `load_harmonised_vcf` runs once (`service.py:979`) and every exporter works from that frame; `study_design` is inferred once and passed in.

**Exclusion reporting is the best in this codebase.** Every per-target loss is reported as its own QC step with a machine-readable reason and a human description of the policy that caused it (`service.py:1150-1212`), including a per-policy explanation of what each duplicate-resolution mode does. `complete_rows` raises rather than writing an empty file when nothing survives.

**`contracts.py` is a practice worth copying.** Each format carries a `FormatContract` recording the tool's real name, what its frequency column means, and a citation to the tool's own documentation. No other module in this review records *why* its output contract is what it is.

**Atomic, collision-checked writes.** `write_dataframe_table` stages to a temporary file and `replace()`s, with a separate uncompressed staging step for gzip targets; output destinations are validated for collisions *before* the VCF is extracted (`service.py:183`, and three tests covering it).

---

## Documentation

`docs/wiki/modules/formatting.md` is accurate, and its Limitations section is honest about the one thing a reader most needs to know: "When LDSC or MiXeR requires study design, it is inferred from count values rather than a separate VCF trait declaration." Three additions worth making:

- State the consequence of that inference: a case-control study whose harmonised VCF carries no case counts will be formatted as quantitative, and LDSC will report observed-scale heritability. That is I2's cross-module risk, and it is the failure a user is least likely to spot.
- Say which targets the `minimum_p_value` clamp can actually affect (MAGMA and GCTA; LDSC uses Z).
- Record why MAGMA and GCTA receive total N while the fine-mapping and MiXeR targets receive effective N (O2).

---

## Suggested order of work

1. **I1** — memoise the identifier resolution on `(identifier_type, duplicate_policy)`. Measured 39.8 s of 48.4 s avoidable on a 5M-row frame with the shipped defaults, and it is the module's dominant cost.
2. **I2** — tighten the binary/quantitative rule beyond "at least one non-null case count", and consider an optional trait-type assertion that is *checked against* the inference rather than replacing it.
3. **Documentation** — the three additions above, especially the LDSC observed-scale consequence.
4. **O1** — decide whether empty PRED-LD partition files should be written and handed to the imputation module.
5. **O2–O5**, and split the two large `resume.py` functions if that file is touched for another reason.

Nothing here changes a value the module writes today. Step 1 is pure cost; step 2 is about making a correct inference harder to get wrong from upstream.
