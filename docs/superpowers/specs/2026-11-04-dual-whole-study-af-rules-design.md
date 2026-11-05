# Dual whole-study AF rules: approved design

The user approved implementation and a local commit after reviewing this design
in chat on 4 November 2026. This supersedes the decision logic in the earlier
whole-study AF acceptance design, not unrelated external-source classifiers.

## Scope and inputs

Keep the existing study-column entry screen: strictly more than 95% of usable
supplied AF values at or below 0.5 triggers whole-study AF preflight. The
non-MAF-like and no-selected-study-column paths retain their existing behavior.
Use the original selected study AF column and the configured independent default
reference, aligned to the original listed study effect allele. Pool the same
eligible, finite, uniquely oriented, non-palindromic biallelic diagnostic pairs
across chromosomes. No new reference file is required. Preserve source rows,
values, row identity and the prepared matches used by chromosome workers.

Let x be unchanged supplied AF, r the reference frequency of the listed effect
allele, and m = min(r, 1-r). The 0.5 boundary is the biallelic minor/major
definition; it is not a new configurable scientific threshold.

## Shared evidence gate

Require at least 1,000 usable pairs and usable coverage >= 0.50. Coverage is
pairs divided by all study rows entering AF preflight after upstream preparation,
before diagnostic exclusions. Equality passes both gates. On failure, report
inconclusive, record the measured values and stop the affected dataset before
chromosome analysis. Do not start an investigation, retry or choose replacement.
Other datasets follow the existing batch failure policy.

## Both complete rules are evaluated

After the shared gate passes, evaluate both complete rules automatically from
the same pooled pairs. Do not short-circuit on rule 1 passing. No user
intervention, extra reference read or second whole-file processing pass is needed.

| Requirement | Rule 1: retain original AF | Rule 2: select replacement |
|---|---|---|
| Correlation | corr(x, r) >= 0.95 | corr(x, m) >= 0.95 |
| Mean absolute difference | mean(abs(x-r)) <= 0.10 | mean(abs(x-m)) <= 0.10 |
| Reference profile | proportion(r <= 0.5) > 0.90 | proportion(r > 0.5) > 0.20 |

Use the configured correlation method (Pearson by default); undefined,
non-finite or constant-vector correlations fail the relevant rule. Mean absolute
difference is not squared. Equality passes correlation and error gates; equality
fails the strict profile gates. A complete rule passes only if all its checks
pass. Rule 2 must compare raw x with m, not min(x,1-x) with m.

Remove the 0.02 directional error-advantage requirement from this whole-study
flow. Do not introduce another error-gap threshold or an informative-subset
classifier. Keep the ordinary external-source/chromosome classifier's existing
behavior and legacy policies; do not silently raise its correlation threshold.

| Rule 1 | Rule 2 | Decision |
|---|---|---|
| pass | fail | retain original AF |
| fail | pass | select explicitly supplied replacement |
| pass | pass | inconclusive; stop affected dataset |
| fail | fail | inconclusive; stop affected dataset |

The both-pass guard must exist and be tested even though the default profile
thresholds are mutually exclusive. Reference-major proportion 0.10 through 0.20
inclusive satisfies neither default profile. The earlier extra exception for
passing EAF agreement metrics while rule 2 passes was explicitly superseded:
only the complete-rule truth table decides. Do not restore that exception.

## Replacement and error categories

When only rule 2 passes, require both resolved --eaf-replacement-file and
--eaf-replacement-column (CLI or canonical YAML). Use shared required-argument
validation, report all missing values together and stop immediately if either
is missing. This is a configuration failure with supported MAF classification,
not scientific inconclusive. Preserve existing early validation of paired options
and invalid input/resource handling. Do not automatically fall back to another
source or continue using rejected original AF. Preserve original MAF provenance,
reference attachment/orientation and existing replacement input validation.

## Configuration, reporting and efficiency

Canonical schema-validated YAML owns all defaults. Add a whole-study correlation
minimum of 0.95 and a strict reference-major minimum of 0.20; reuse the shared
whole-study count, coverage, mean-error and minor-profile controls. A separate
whole-study correlation key is necessary because the existing correlation key
also controls unchanged ordinary external-source classification at 0.90.
Register new keys in effective configuration, exports, preflight identity and
logging. Reject out-of-range, non-finite, boolean and otherwise invalid values
through existing policy validation; no module hardcoded defaults.

Record both full rule results, each failed check, counts, coverage, both raw
correlations, both errors, both profile proportions, thresholds and final action
in canonical evidence and study reports, including failures. Existing consumers
must receive consistent status/source/provenance. Preserve shared progress,
screen recording, required-argument and checkpoint infrastructure. Reuse existing
reference preparation and vectorized aggregation. No protected adapter edits,
new dependencies or broad restructuring.

## Scientific interpretation and limits

These are user-selected QC heuristics, not proof of regression coding or of every
row's frequency semantics. EAF is frequency of the listed effect allele and may
be below or above 0.5. Population differences near 0.5 can make EAF and MAF fits
similar. Record this known limitation without adding a new decision condition.
No 100–200-million-row benchmark or production GWAS rerun is implied by unit
tests. Validate behavior with small representative real interfaces, exact
boundaries, missing/invalid inputs and synthetic controls before committing.

Supporting definitions and API verification:

- [GWAS-SSF definitions](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.1.0.tex)
- [Winkler et al., Figure 4 and AF QC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083217/)
- [Polars correlation API](https://docs.pola.rs/py-polars/html/reference/expressions/api/polars.corr.html)
