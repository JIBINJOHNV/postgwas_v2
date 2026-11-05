# Optional EAF replacement source

Use the paired options below when a study supplies an internal frequency column
that may be minor-allele frequency (MAF):

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --resource-directory resources \
  --output-directory results \
  --eaf-replacement-file 'panel_{build}_{chromosome}.tsv.gz' \
  --eaf-replacement-column EUR
```

The equivalent full-run YAML is:

```yaml
modules:
  harmonisation:
    eaf_replacement_file: 'panel_{build}_{chromosome}.tsv.gz'
    eaf_replacement_column: EUR
```

Both canonical values default to `null`. Set both together; a missing partner
reports its public option and canonical YAML path. Explicit command-line values
override YAML. Export these settings through the normal configuration command:

```console
postgwas config export --module harmonisation --style full --output harmonisation.yaml
```

The file may contain all chromosomes or use `{build}` and `{chromosome}` path
placeholders. Relative CLI and run-YAML paths resolve from the current working
directory, matching existing run resource settings; sample-sheet row paths
continue to resolve from the sample-sheet directory. Other placeholders, format specifiers and filename prefixes are
invalid. The replacement value column must contain frequencies for the allele
declared by `external_eaf_mapping.effect_allele` (ALT in the packaged mapping).
The structural chromosome, position and allele columns and delimiter come from
`external_eaf_mapping`; no sample-sheet column is added.

PostGWAS first tests the selected internal frequency against `default_eaf`.
Accepted EAF stays in use. Replacement is provisionally selected only when the existing
comparison panel supports the MAF interpretation. A confirmed MAF input without
an explicit replacement stops that dataset. Insufficient or conflicting evidence
continues to stop the dataset. The replacement table is then checked by the
separate whole-study gate below before its EAF can be adopted. It does not
supply the initial MAF decision or palindrome orientation. An ordinary sample-sheet external EAF source keeps its existing
source-selection and comparison rules; frequency attachment uses the same
validated allele-frequency join described below.

For a MAF-like study column, require at least 1,000 pooled usable pairs and
coverage >= 50% of all rows entering AF preflight after upstream preparation.
Diagnostic exclusions remain in that fixed denominator. Count and coverage
equality pass; failure is inconclusive and stops before either rule or replacement.

With `x` unchanged supplied AF, `r` reference AF aligned to the original listed
study effect allele, and `m = min(r,1-r)`, evaluate both complete rules:

| Requirement | Rule 1: retain original AF | Rule 2: select replacement |
|---|---|---|
| Correlation | `corr(x,r) >= 0.95` | `corr(x,m) >= 0.95` |
| Mean absolute difference | `mean(abs(x-r)) <= 0.10` | `mean(abs(x-m)) <= 0.10` |
| Reference profile | `proportion(r <= 0.5) > 0.90` | `proportion(r > 0.5) > 0.20` |

Rule 2 never folds the supplied AF. Both rules are measured even if rule 1
passes. Only rule 1 passing retains original AF; only rule 2 passing provisionally selects
the explicit replacement. Both passing or both failing is inconclusive and
stops the dataset. A partial EAF fit with a full rule-2 pass selects replacement.
No 0.02 directional error-advantage margin applies to either whole-study rule.
Reference-major fractions 0.10 through 0.20 inclusive meet neither default profile.

All whole-study gates are configured under `modules.harmonisation.policies.eaf`:

```yaml
modules:
  harmonisation:
    policies:
      eaf:
        maf_reference_min_overlap: 1000
        maf_reference_min_overlap_fraction: 0.50
        study_reference_min_correlation: 0.95
        maf_reference_max_mean_absolute_difference: 0.10
        study_reference_min_minor_fraction: 0.90
        study_reference_min_major_fraction: 0.20
```

Use finite proportions from 0 to 1, not `90`. Correlation and error equality
pass; profile equality fails. AF exactly 0.5 counts toward the minor fraction.
Undefined, non-finite or constant-vector correlations fail the relevant rule.
The shared method remains Pearson by default, with Spearman supported.
The separate ordinary external classifier retains
`eaf.maf_reference_min_correlation: 0.90`,
`eaf.maf_reference_error_margin: 0.02` and its upper directional
`eaf.reference_minor_fraction_cutoff: 0.95`.

These are user-selected QC heuristics. A deterministic 10,000-variant MAF
control with 94% minor reference effect alleles passes the gates (correlation
0.9523683946; minor fraction 0.94; direct mean absolute difference 0.012).
Passing does not prove the semantic provenance of the supplied column. Population
differences near 0.5 can make EAF and MAF fits similar; this known limitation
does not add a decision condition. See the
[method and configuration discussion](configuration.md) for the fixed-denominator
contract and this limitation.

Messages distinguish study EAF retention, provisional replacement selection,
replacement approval and a stopped dataset. Supported MAF without a replacement reports both required options and
their YAML keys. It stops that dataset before chromosome analysis; other
datasets continue under the existing batch failure policy.

## Whole-study replacement validation

After replacement resources and chromosome projections are staged at dataset
step 08, PostGWAS uses the existing allele-frequency join and validated
default-reference orientation to measure the candidate once across the study.
Let A be unchanged original study AF, B be replacement frequency expressed for
the **original listed study effect allele**, and C = min(B, 1-B). Every gate
must pass before chromosome workers can adopt the source:

| Measurement | Packaged requirement | Existing policy under `policies.eaf` |
|---|---|---|
| Usable pairs n | >= 1,000 | `maf_reference_min_overlap` |
| Coverage n/N | >= 0.50 | `maf_reference_min_overlap_fraction` |
| corr(A,C) | >= 0.95 | `study_reference_min_correlation` |
| mean(abs(A-C)) | <= 0.10 | `maf_reference_max_mean_absolute_difference` |
| Fraction B > 0.5 | > 0.20 | `study_reference_min_major_fraction` |

N remains all study rows entering original AF preflight after upstream
preparation. Count, correlation, error and major fraction use the same finite,
valid, uniquely aligned, non-palindromic biallelic pairs. Diagnostic exclusions
stay in N and do not remove study/output rows. Count and coverage failure stops
before correlation or ranking; those uncomputed metrics are reported as
unavailable. Equality passes count, coverage, correlation and error; equality
fails the major fraction, and B exactly 0.5 is not major. Original A is never
folded, changed or used to infer strand. The later canonical ALT direction
does not redefine the original-EA major fraction. The worker's separate
`external_min_match_fraction` coverage requirement still applies per chromosome.

The shared `maf_reference_correlation_method` is Pearson by default. Pearson
combines centered sufficient statistics; Spearman uses exact global ranks with
ties, never chromosome-local ranks or an average of chromosome correlations.
No extra directional-error margin, reference panel, retry or fallback is added.
Any failed numeric gate records every failed requirement as `inconclusive`
and stops the dataset before workers. A changed source/work identity after
measurement records `invalid` with `replacement_inputs_changed`; earlier
I/O/schema failures retain their original exception and measured original
evidence without fabricating replacement metrics.

`eaf_reference_validation` preserves the original dual-rule decision.
`eaf_replacement_validation` records the separate status, source file/column/
mapping and identities, fixed N, paired input rows, usable count/coverage,
diagnostic exclusions, method, all thresholds, metrics and failure reasons.
Only final `passed` evidence with its row/source/policy-bound approval contract
allows adoption. A numeric helper's earlier passing log does not override a
later identity failure. The terminal shows a later replacement summary after
the initial provisional decision; canonical logs and HTML retain both.
Same-source file aliases reuse matching prepared evidence and explicitly
record `same_source_evidence_reused` and `reused_default_chromosomes`; this is
reused evidence, not independent confirmation. A different path alone cannot
establish scientific independence.

Passing supports consistency under the declared counted-allele mapping; it
does not prove per-row allele direction or recover measured study EAF.
For example, repeating A=[.1,.2,.3,.4] and B=[.1,.8,.3,.6] gives correlation 1,
MAE 0 and major fraction .5. Replacing B by 1-B gives exactly those same
metrics. This mapping limitation introduces no additional decision gate.
The [GWAS-SSF contract](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.1.0.tex)
still requires the frequency to refer to the declared effect allele;
[frequency QC guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083217/)
does not calibrate these user-selected numerical thresholds.

The selected frequencies are a population-reference proxy, not recovered
study-measured EAF. Choose the intended population and source genome build.
The [GWAS-SSF field definitions](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex)
define frequency relative to the listed effect allele; the effect allele need
not be minor. [1000 Genomes](https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/)
reports ALT frequencies, including separate values for multiple ALT alleles.
The [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf) defines
one AF value per ALT allele. Diagnostic MAF folding remains restricted to
biallelic reference sites and does not change when replacement is selected.

Direct matches retain the listed frequency. A validated single allele pair
retains the existing swapped `1-AF` calculation, including supported sequence
indels. For a complete SNP site with one common REF and distinct A/C/G/T ALT
alleles, the shared join can calculate the REF frequency once as
`1 - sum(all distinct ALT frequencies)` and use it for swapped matches. For
example, `REF=A`, `ALT=G,T`, and ALT frequencies `0.30,0.10` give REF frequency
`0.60`; subtracting only `0.30` would give the incorrect `0.70`.

Complete site membership in the supplied reference remains a declared input
assumption. All components must describe the same allele set and population.
Each required frequency must be present, finite and within 0 to 1, and the
sum must leave a valid REF frequency. Exact duplicates follow the configured
duplicate policy and count once when retained. A discarded or conflicting
component cannot silently disappear from the sum. Different REF alleles,
unsupported alleles, and multiple overlapping indel records do not establish
an eligible SNP allele set merely by sharing a chromosome and position.

An unresolved swapped assignment receives a missing generated EAF; original
study MAF and reference values remain unchanged. It stays in the denominator
of `eaf.external_min_match_fraction` and follows the existing required-field
missing policy. With the packaged `final_check.on_missing: reject`, it is
removed and recorded rather than assigned a guessed frequency. The existing
coverage, rejection and explicit duplicate-failure policies still apply.
Reference I/O or configuration failures remain errors.

This calculation applies to reference-frequency attachment. The bundled
GWAS2VCF adapter remains unchanged and still uses `1-AF` when it reverses
alleles. At a multiallelic site, that complement can include other alleles
rather than the required opposite allele. No additional export check is made,
so complete-site attachment does not provide unrestricted multiallelic VCF
allele-swap support. The existing `vcf.liftover_swap` policy remains unchanged.

Site calculation reuses the loaded reference and existing study join. It does
not rerun the study-wide MAF/EAF decision or read the reference again for the
calculation. The aggregation adds CPU and temporary-memory work. Runtime and
memory impact require measurement on the intended data and environment. This
calculation does not create or rewrite reference files.

A supplied replacement specification must match at least one existing nonempty
file, even if a particular dataset later retains its internal EAF. Exact active
build/chromosome coverage and table contents are validated during dataset
preparation. A shared whole-genome replacement uses the same bounded staging
and worker memory accounting as existing external EAF tables.

Replacement diagnostics process one narrow chromosome and its reference/merge
working set at a time, then replay temporary paired Parquet files; they do not
repeat reference parsing/merging for numeric passes. Distinct candidates are
attached again during workers because diagnostics omit ineligible output rows.
The premerge preparation step has no configured RAM-admission check: a large
chromosome or reference can exceed memory despite this structural bound.
Exact global Spearman materializes eligible pairs/ranks only after its existing
execution-budget estimate; a null direct-library budget is unbounded. That
planning check is not an OS memory limit or a disk-spill ranking implementation.
Temporary diagnostic files are removed on success/failure. Full-workflow
100–200-million-row capacity remains unverified; numeric-only and smaller shared
merge timings cannot establish reference I/O, staging or full-run capacity.

Before a completed direct command can resume, its replacement template is
expanded and the resulting exact files join the durable checkpoint input
fingerprints. This applies to files outside `resources.root`, YAML settings,
`--option VALUE` and `--option=VALUE`. A changed tracked source invalidates the
completed checkpoint under the canonical restart policy. Original study MAF,
replacement source and final EAF interpretation remain distinct in provenance.
