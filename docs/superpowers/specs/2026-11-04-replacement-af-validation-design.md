# Whole-study replacement AF validation: implemented design

Implemented under the user's subsequent implementation authorization. This
document specifies the post-selection gate added to the original-study rules
implemented in `807f57b`. Initial inspection baseline: `15fd11c`; unrelated
exponentiated-SE work was subsequently committed as `e943add` and preserved.
Runtime implementation and independent review are complete. Final scoped
verification passed: 782 workspace cases and 317 clean-candidate cases. The
[implementation plan](../plans/2026-11-04-replacement-af-validation.md#final-implementation-and-validation-record-2026-11-04)
records exact commands, test/file scopes, all six ledger rulings, baseline
failures and the authorized local-only commit procedure; no push is implied.

## Policy and scope

Keep the existing original-study entry screen, default-reference comparison and
two-rule decision table. After sufficient original evidence, evaluate both rules:
only rule 1 passes -> retain original AF; only rule 2 passes -> provisionally
select the explicit replacement; both pass or both fail -> inconclusive and stop.
Missing replacement file/column is an actionable configuration error. A supplied
but invalid file/column is an invalid argument, not a missing argument.

Only provisionally selected replacements enter the new gate. Use the existing
external-frequency merge and its configured allele mappings, normalization,
duplicate handling, direct/swapped alignment and complete-site multiallelic
calculation. Do not implement another attachment algorithm. Identical default
sources reuse validated aligned frequencies and evidence where the complete
comparison contract matches.

Let A be the unchanged original study AF, supported as MAF by the first decision;
B be replacement frequency expressed for the original listed study effect
allele; C = min(B, 1-B). The new whole-study decision requires all five gates:

| Measurement | Requirement | Existing canonical policy under `modules.harmonisation.policies.eaf` |
|---|---|---|
| Usable comparison count n | n >= 1,000 | `maf_reference_min_overlap` |
| Coverage n/N | >= 0.50 | `maf_reference_min_overlap_fraction` |
| Correlation of A with C | >= 0.95 | `study_reference_min_correlation` |
| Mean absolute difference of A and C | <= 0.10 | `maf_reference_max_mean_absolute_difference` |
| Fraction of B > 0.5 | > 0.20 | `study_reference_min_major_fraction` |

Reuse `maf_reference_correlation_method`: Pearson is the canonical default;
Spearman must retain global rank semantics when selected. Do not silently switch
methods or average chromosome correlations. No directional error-advantage
margin or extra reference panel is added. Reuse these existing policies rather
than introduce duplicate defaults; update their YAML help to cover both uses.

N is the fixed number of study rows entering AF preflight after upstream
preparation, before diagnostic exclusions. All three frequency metrics and n use
the same eligible, finite [0,1], unique, non-palindromic biallelic paired rows.
Diagnostic exclusions do not discard source rows and remain in N. Count/coverage
failure stops before correlation or ranking. Equality passes count, coverage,
correlation and MAE; exactly 20% major fails; AF exactly 0.5 is not major.

Multiallelic attachment stays supported by the existing merge. General
multiallelic MAF is not `min(B,1-B)`; these rows keep the existing exclusion from
MAF diagnostics. Do not add an output-row exclusion because of this diagnostic.

If any numeric evidence gate fails, record replacement validation as
`inconclusive` and stop the dataset before scientific chromosome workers start.
Report each failed requirement. Do not retry, fall back to original AF as output
EAF, select another source or loosen thresholds. Preserve original AF provenance.
I/O and malformed/scientifically
incompatible input failures retain their existing error classifications.

## Orchestration and provenance

Place the new gate after shared replacement staging during dataset preparation,
before preparation succeeds or workers launch. Source selection is provisional
until the gate passes. Keep original `eaf_reference_validation` evidence intact;
add distinct `eaf_replacement_validation` evidence and require its `passed` status
at replacement adoption. Failure is also reported through the existing canonical
study-failure and progress infrastructure, with no completed scientific output.

The attachment merge expects strand-aligned keys. Prepare narrow temporary keys
from validated default-reference orientation evidence; do not feed raw reverse-
complement keys to its direct/swap-only join. Never use original MAF or candidate
replacement AF to establish palindrome/strand direction. Retain row identity and
the original study alleles. Restore B to original study-EA direction for the
global major fraction; the later canonical ALT direction must not change that
fraction. Do not run effect conversions twice or modify BETA, Z or original AF.

Read staged chromosome projections under existing resource bounds, reuse matching
and complete-site calculations, and retain/reuse their results where safe. Pool
Pearson sufficient statistics with numerically stable chunk-combination and one
fixed denominator. Exact global Spearman requires global ranking under existing
resource controls; chunk-local ranking is invalid. Do not claim bounded memory
for the whole pipeline or 100-200-million-row capacity from unit tests.

Record source file/column/mapping, original-MAF evidence, counts/exclusions, N,
method, thresholds, correlation, MAE, major count/fraction, outcome and failure
reasons in canonical logs, terminal and HTML. Bind reused evidence to immutable
source identities, allele mappings, policies, software and row identity. Changed
inputs/configuration must invalidate the checkpoint through the existing policy.
Retain declared/formatted source and staged pathnames separately from resolved
identities, re-resolve the same pathnames before final approval and adoption,
and pin the verified physical attachment path across intervening strand work.
A retargeted alias must not supply new frequencies under old evidence. This
path pinning is not a filesystem lock against writes to the physical file.

## Interpretation and supporting sources

Passing means replacement frequencies are consistent with the original study
MAF under the declared allele mapping and meet the configured major-frequency
profile. It does not prove per-row allele direction or recover measured study
EAF: folding B and 1-B gives the same MAF. A mislabeled MAF source can acquire
values above 0.5 through a mistaken complement. Keep the source's declared
counted-allele contract and existing independent orientation requirements.
Identical replacement/default sources supply reused, not independent, evidence.
The numerical thresholds are user-selected QC policies, not calibrated error rates.

- [GWAS-SSF effect-allele and EAF definitions](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.1.0.tex)
- [VCF AF definitions, one value per ALT](https://samtools.github.io/hts-specs/VCFv4.5.pdf)
- [Winkler et al., study-level frequency QC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083217/)
- [Polars correlation methods](https://docs.pola.rs/py-polars/html/reference/expressions/api/polars.corr.html)

## Implementation rulings and verification limits (2026-11-04)

- The active feature checkout is reused, with task snapshots isolating unrelated
  user changes; one final local commit is authorized, with no push implied.
- Numeric `passed` becomes approval only after final input checks issue the
  row/source/policy/software-bound contract. Numeric failure is `inconclusive`
  without a contract. A late identity change is `invalid`, preserving metrics
  with `replacement_inputs_changed`; earlier I/O/schema exceptions keep their
  original classification and original AF evidence.
- Same-source aliases/columns/mappings record reused evidence. Distinct file
  paths do not establish independent scientific provenance.
- Preparation holds one narrow chromosome plus its reference/merge working set
  and replays temporary paired Parquet. It does not enforce configured RAM
  admission before that merge. The later worker scheduler estimates different
  work and is not a preparation-capacity guarantee. Shared preparation admission
  would require a separately reviewed scheduling extension.
- Exact global Spearman retains the existing execution-budget planning check.
  Null direct-library budgets remain unbounded; neither the estimate nor the
  per-chromosome structure is an OS memory limit, and no spill ranking exists.
- The A=[.1,.2,.3,.4], B=[.1,.8,.3,.6] characterization gives corr=1, MAE=0,
  major fraction=.5 for both B and 1-B. It protects the declared limitation,
  with no extra gate. An unequal 10/990-pair adversary has each local Pearson
  correlation 1 but pooled correlation 1/199; global pooling must fail .95.

Final verification recorded 782 passing tests with 78 existing warnings in
277.14 s across 25 workspace owners, and 317 passing tests in 108.58 s across the
nine affected owners in the exact clean candidate. The broad run preceded only
the final path-constant restoration for an existing local test importer; the
clean-candidate run includes that restoration. No runtime source or test function
body changed for that export. Candidate tests generate their own synthetic data
and do not read the local manifest or private GWAS payload. Independent review
approved the alias-retarget correction, fixture portability and restored export.

Final static checks passed: 15 Python owners compiled; the three reviewed
numerical-function ASTs and six canonical AF defaults were unchanged; no protected
adapter edits occurred. YAML changes are help/origin text only, with no new key,
default, flag or cutoff. Canonical policy synchronization and the 50-page wiki
check passed. Environment: Python 3.10.21, Polars 0.20.31, NumPy 1.23.5,
bcftools/htslib 1.23.1 on Darwin arm64 with 12 CPUs and 64 GiB RAM.

Unrestricted pytest still has exactly the same 37 baseline collection errors
(9.89 s), with no introduced or removed error. The plan preserves all 37 paths.
The separate documentation CLI check remains red due to a pre-existing command-
extraction false positive in valid installer `--name postgwas --jobs 2` examples;
the installer supports `--jobs`, and validator repair is outside this task.
The broad workspace suite includes a pre-existing tracked exponentiated-SE test
that imports an ignored concordance helper, so no fresh-clone claim is made for
the entire tracked suite. The clean candidate verifies the nine affected owners.
All six ledger decisions and their costs, the 21 exact owned files, commands,
review outcomes and local-commit procedure are permanently recorded in the plan.

Numeric-only generated Pearson pairs at 100M/200M took 4.5046/7.5457 s with
20-ms sampled process RSS peaks 207,536,128/209,420,288 bytes. They performed no
reference I/O, allele join, global ranking or full workflow. A separate shared
merge benchmark used 100K/1M study rows and 110K/1.1M reference rows with direct,
biallelic-swap and complete multiallelic REF sums; expected AF matched at 1e-15.
It took 0.2182/1.9767 s and sampled 323,158,016/1,577,566,208 bytes RSS.
These layer-specific measurements do not verify full-workflow 100–200M capacity,
OS peak memory, production I/O/spill behavior or source-panel completeness.
