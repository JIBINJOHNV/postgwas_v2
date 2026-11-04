# Whole-study AF acceptance thresholds — proposed design

Status: implemented with the directional error margin unchanged at `0.02`. Scoped validation and its limitations are recorded in the implementation plan.

## Requested behavior and scope

For the whole-file validation of a selected, MAF-like study AF column, retain the original AF only when all of these conditions pass:

| Condition | Planned default and operator |
|---|---|
| Same-effect-allele frequency correlation | `correlation >= 0.90` |
| Reference fraction with aligned effect-allele AF at or below `0.5` | `fraction > 0.90` |
| Usable diagnostic pairs | At least `0.50` of study rows entering preflight, and at least the existing absolute minimum of `1,000` |
| Mean absolute allele-aligned AF difference | Existing maximum, currently `0.10` |

The interpretation of “below 0.5” retains the existing inclusive `AF <= 0.5` metric. Exactly `0.5` counts in the numerator. This is an explicit boundary choice for review; a literal `< 0.5` definition would require a different metric and boundary tests. Correlation and proportions use values from `0` to `1`; the configuration value is `0.90`, not `90`.

Pool the existing eligible matched variants across all chromosomes. Calculate the correlation and fraction from that pooled set, rather than averaging chromosome correlations or fractions. Use the same valid, uniquely oriented, non-palindromic, biallelic comparison pairs already prepared by the current implementation.

For the new overlap fraction, the numerator is the pooled usable diagnostic-pair count and the denominator is all study rows entering whole-file AF preflight after existing upstream preparation. The denominator is fixed before diagnostic exclusions: palindromic, ambiguous, non-biallelic, missing-frequency and invalid-frequency rows cannot contribute pairs, but remain in that denominator. This is usable AF-comparison coverage, not merely coordinate/allele-key coverage. Exactly `0.50` passes. Retain the existing `1,000` absolute floor, so a smaller dataset must still meet that floor even if every row matches.

The new overlap-fraction and reference-minor-fraction gates apply to the current whole-study native-AF acceptance path. The user explicitly requests the existing shared correlation default to change to `0.90`, affecting conditional MAF inference and ordinary external-source classification. The latest instruction excludes the directional error-margin change, so its existing `0.02` default is preserved. Ordinary external AF attachment and classification are not moved to the whole-file stage. The WLM external-source error-margin failure is not resolved by this plan: its unchanged required gap is `0.02`, and its recorded chromosome-22 gap of approximately `0.018279` remains below that margin.

## Approaches

1. Reuse the existing shared correlation setting, as now explicitly requested, and add a separate whole-file overlap fraction. Recommended: the existing integer overlap policy keeps its count semantics and the new fraction has a clear denominator. Keep one correlation setting rather than creating a duplicate with the same new default.
2. Change the existing integer overlap setting to accept proportions or continue adding a direct-only correlation setting. Rejected: either changes the units of an established option or duplicates the explicitly requested shared correlation setting.

This revision supersedes the earlier proposed `eaf.study_reference_min_correlation` key, which was never implemented. The separate reference-minor-fraction acceptance setting is still needed because the existing directional cutoff has the opposite purpose.

## Configuration

Define these in the canonical policy registry, `src/postgwas/config/defaults/modules/harmonisation.yaml`:

- Add `eaf.maf_reference_min_overlap_fraction`, default `0.50`; inclusive whole-file usable-overlap acceptance.
- `eaf.study_reference_min_minor_fraction`, default `0.90`; strict acceptance.
- Change `eaf.maf_reference_min_correlation` from `0.80` to `0.90`; shared inclusive correlation requirement.
- Preserve `eaf.maf_reference_error_margin: 0.02`; no directional-margin behavior change is authorized.

The two new proportions use the existing finite numeric `[0, 1]` validator, documented help, provenance, EAF ordering, and common export list. Preserve the absolute `maf_reference_min_overlap` default of `1,000` and the mean-difference maximum of `0.10`. This YAML supplies both validation definitions and default values. Do not duplicate defaults in Python or create a second defaults file. Full and common configuration exports must expose and reload the new and changed values.

The shared `HarmonisationConfig.policies` model boundary also invokes the existing registry validator, so both full-run and module-only public configuration validation reject invalid overrides before analysis. It preserves the original override mapping and translates registry errors through the shared configuration-error formatter; it adds no per-key validation, duplicate defaults or optional analysis dependency.

```yaml
modules:
  harmonisation:
    policies:
      eaf:
        maf_reference_min_overlap: 1000
        maf_reference_min_overlap_fraction: 0.50
        maf_reference_min_correlation: 0.90
        maf_reference_max_mean_absolute_difference: 0.10
        maf_reference_error_margin: 0.02
        study_reference_min_minor_fraction: 0.90
```

Keep `eaf.reference_minor_fraction_cutoff: 0.95` for its existing upper-limit directional-evidence purpose. Keep `eaf.external_min_match_fraction: 0.80` for chromosome-level external attachment coverage; it is not the requested diagnostic overlap fraction. The existing `0.02` margin does not participate in a directly accepted whole-study EAF decision.

## Decision, provenance, and reporting

- Extend `_compare_frequency_pairs` rather than adding another matching or classifier implementation. Add an optional keyword `study_variant_count=None`; require a positive total for the whole-study `accept_direct_eaf=True` path and pass `df.height` once from `validate_study_allele_frequency`. Ordinary external calls keep their existing absolute overlap requirement without a new per-chromosome fraction gate.
- Preserve `reference_effect_allele_minor_fraction`; calculate no new reference join or reread.
- Record `study_variant_count`, `observed_overlap_fraction`, `minimum_overlap_fraction`, `minimum_overlap` (the absolute floor), and `minimum_reference_effect_allele_minor_fraction`. Existing direct/folded correlation evidence can both use the shared `0.90` threshold; no separate direct-correlation policy is needed.
- Require both the absolute overlap floor and `paired.height / study_variant_count >= configured_fraction` before using the pairs for either whole-file EAF acceptance or MAF inference. Insufficient fractional overlap returns inconclusive with `insufficient_overlap_fraction` and cannot authorize replacement. Compare the fraction directly rather than using a floating-point `ceil(fraction * count)` that can reject exact boundary matches under overrides such as `0.14` of `100`.
- A new fraction-gate failure on otherwise passing direct agreement must return inconclusive with `reference_minor_fraction_not_above_minimum`. It must not alone trigger MAF replacement or fall through to a weaker EAF acceptance.
- Existing replacement remains possible only when direct agreement itself fails and the separate MAF classifier supports MAF; a replacement file and column remain required. Other unresolved cases stop before chromosome analysis.
- Display observed correlation, observed reference-minor fraction, thresholds, operators, usable-pair counts and failed condition in existing shared reporting. Preserve original AF values and all source rows during preflight.
- Reuse the dataset decision in chromosome workers. Preserve ordinary row QC, orientation, and explicitly selected replacement behavior. Ordinary external/folded-MAF comparisons intentionally inherit the changed shared correlation default and unchanged error margin, while their absolute overlap behavior remains unchanged.
- Verify that resolved policy changes participate in existing checkpoint identity and invalidate old decisions under the configured restart policy. Do not introduce a separate checkpoint system.
- Guard undefined correlations explicitly when either paired vector is constant, for both direct and folded comparisons. Reuse the min/max guard pattern in `src/postgwas/modules/harmonisation/concordance/analysis.py`; do not add a guessed numerical epsilon or change valid nonconstant-vector decisions.

## Scientific interpretation and limitation

These are user-selected QC heuristics, not a validated proof that a column is EAF. A deterministic synthetic control with 10,000 variants, 94% reference minor effect alleles, and a supplied column equal to MAF gives correlation `0.9523683946`, reference-minor fraction `0.94`, and direct mean absolute difference `0.012`. It passes the proposed gates despite containing MAF irrespective of effect-allele identity. Document this limitation; do not describe the change as resolving all EAF/MAF ambiguity or as scientifically calibrated.

The [GWAS-SSF definition](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.1.0.tex) ties EAF to the listed effect allele. It does not validate these numerical thresholds. The control above uses the formulas in the current `_compare_frequency_pairs` implementation; no additional study collection is needed for it.

## Validation and constraints

Cover inclusive correlation, inclusive overlap and strict reference-minor-fraction boundaries; odd row totals; the absolute floor; fixed denominator across chromosomes and diagnostic exclusions; missing/invalid values; constant vectors; no overlap; preserved study values; one pooled decision; reporting; configuration export/override; checkpoint identity; and the explicitly changed external-source thresholds. Boundary examples: `500/1000` passes the fraction gate, `499/1000` fails it, and `14/100` passes a `0.14` override without rounding the count upward; the absolute floor is tested separately. Reuse the existing real-record study-EAF fixture and small synthetic controls. No 100–200-million-row scalability claim is made: the existing whole-study preparation is not a newly bounded-memory implementation.

Protected adapters and unrelated files remain outside scope. Implementation requires the repository change-cycle checks, focused scientific regression tests, complete cumulative diff review and a report of any unresolved limitations.

## Baseline verification finding

The initial configuration-export, study-reporting and AF-failure-reporting baseline produced `54 passed, 1 failed` before runtime changes. `test_inconclusive_af_failure_retains_observed_comparison_for_reporting` expected `correlation_unavailable` for constant study AF; Polars `0.20.31` supplied a finite artifact (`-1.2593779783375902e-07`), so the old comparator reported `correlation_below_minimum`. Both outcomes stopped the dataset, but the diagnostic was incorrect. The implemented equal-min/max guard now reports unavailable direct and folded correlations correctly; the existing scientific failure expectation was retained.
