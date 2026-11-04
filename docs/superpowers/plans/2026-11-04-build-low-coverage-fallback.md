# Genome-build low-coverage fallback implementation

Follow-up decision (2026-11-04): the user requested
`low_coverage_fallback_enabled: true` as the canonical default. This supersedes
the original opt-in/default-disabled requirement below. Explicit `false` still
disables fallback; the scientific thresholds and decision order are unchanged.

The user authorised implementation of the proposed package-wide fallback and
asked whether 10,000 means matches to either build. The required interpretation
is **10,000 distinct SNP loci supporting only the winning build**, not the sum
of matches to either build. This document is the implementation specification.

## Global constraints

- Follow AGENTS.md. Preserve unrelated `_to_delete/` and `tmp/` directories.
- Do not modify protected harmonisation adapters.
- Retain normal automatic acceptance and declared-build behavior. The fallback
  is opt-in and disabled by the canonical YAML default.
- The fallback may waive reference coverage only when it is the sole failed
  existing automatic gate. All other normal gates must pass.
- Require winner matches / coordinate-testable study rows >=0.95, competitor
  matches / the SAME denominator <=0.05, and >=10,000 distinct winner-only SNP
  loci. Coordinate-testable rows include allele-incompatible coordinate hits.
  Preserve existing row-based fractions; explicitly test duplicate handling.
- Unique evidence counts canonical (chromosome,position) once, includes only
  compatible single-nucleotide A/C/G/T substitutions with different alleles,
  and excludes a locus if ANY study evidence there matches the other build.
  Both-compatible matches are not exclusive evidence. Do not filter study rows.
- Veto fallback if any chromosome has >=100 unique competitor-only SNP loci
  AND competitor matches / coordinate-testable study rows on that chromosome
  >=0.95. Unassessed or sparse chromosomes are not certified consistent.
- Keep per-row/strand/reference evidence semantics and downstream interfaces.
  Add structured acceptance-path, distinct-locus, per-chromosome and fallback
  decision evidence to the returned dictionary/manifests, canonical log, and
  study report. Make failure messages explain which fallback gates failed.
- All defaults/thresholds must be canonical schema-validated YAML policies;
  do not add independent parser defaults or hardcode scientific cutoffs.
  Use flat build policy names consistent with the existing registry:
  `low_coverage_fallback_enabled` false,
  `low_coverage_min_match_fraction` 0.95,
  `low_coverage_max_competitor_fraction` 0.05,
  `low_coverage_min_unique_loci` 10000,
  `low_coverage_chromosome_min_competitor_loci` 100,
  `low_coverage_chromosome_min_competitor_fraction` 0.95.
- Fallback inference compares exactly two candidate references. Preserve normal
  behavior with other reference counts; when automatic fallback is enabled,
  reject an incompatible reference count with an actionable error.
- Sources support comparative mapping but not these exact thresholds:
  https://cloufield.github.io/gwaslab/InferBuild/ (10,000-match warning),
  https://doi.org/10.1093/bioadv/vbaf262 (10,000 sampled rows), and
  https://github.com/neurogenomics/MungeSumstats/blob/master/R/get_genome_build.R.
  Document the policy as experimental dominant-build inference. A 95% match
  fraction is not a posterior probability; small mixed components can pass.

### Task 1: Implement, test and document the complete fallback

Read the relevant owning files and complete call paths before editing. Extend
the existing inference join/evidence collection rather than reading references
twice or adding a parallel allele-matching implementation. Reuse the canonical
configuration/schema/export mechanism and existing logging/UI infrastructure.

Write and run failing real-interface tests first. Cover default-disabled normal
behavior, symmetric accepted low-coverage builds, exact and just-failing 95/5
boundaries and 10,000 unique SNP floor, matches to both builds, incompatible
allele coordinate hits, duplicate study/reference rows, same-locus opposing
alleles, non-SNP evidence, chromosome veto boundaries, normal acceptance despite
veto (scope preservation), declared-build bypass, multiple failed normal gates,
zero/missing evidence, valid/invalid policies, canonical export/log/report and
preserved strand evidence. Use small parameterized evidence floors only for
isolated behavior tests; include the actual 10,000 default boundary once.

Update canonical harmonisation YAML required/default entries and registry
lifecycle lists, affected exact exported-key tests, configuration docs, and
study reporting/QC summaries if needed. Provide a complete opt-in YAML example
using existing `--run-config` support. Do not add an unnecessary CLI option.
Ensure the effective policies contribute to existing checkpoint identities.

Expected: focused tests pass, existing build/strand/report/configuration and
checkpoint behavior remain valid, disabled defaults preserve prior scientific
results, no protected/unrelated files change. Report red/green commands and
outcomes, complete cumulative diff review, scientific limitations and checklist.
Do not commit or push; root integrates/reviews the cumulative implementation.

### Task 2: Independent review and real-data validation

Root reviews the cumulative change and runs the updated production function on
the retained independently documented GRCh37/GRCh38 datasets used in the prior
offline experiment. Confirm fallback recovers the predefined ADHD/GRCh38 25%
subsets with enabled configuration, preserves disabled/default results and
rejects the duplicate and concentrated mixed-build controls. Report the known
accepted minor-mixture limitation. Do not alter thresholds to fit outcomes.

Run the relevant test suite and repository test command; record pre-existing
out-of-scope failures without modifying unrelated code. Run applicable lint/
configuration checks, inspect the whole diff, and request independent review.
Any material finding requires a focused regression and re-review. Final handoff
must include activation instructions, checks, limitations and accurate status.
