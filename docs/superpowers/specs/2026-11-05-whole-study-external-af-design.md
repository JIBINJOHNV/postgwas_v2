# Whole-study initial external allele-frequency validation

The user requested correction of the no-mapped-study-AF path after the EGG
recovery log reproduced the old chromosome-level margin failure on current
commit `1ea6305`. This extends the earlier whole-study design to the initially
selected external source. It does not change the approved scientific thresholds.

## Required behavior

When no study AF column is mapped, retain that finding and select the configured
`external_eaf_file`/column. After existing chromosome partitioning and shared
external-file staging, validate this source once across the dataset before any
chromosome worker starts. Store separate `eaf_external_validation` evidence;
do not relabel proxy reference frequencies as measured study AF or as an
explicit replacement. Native-study AF and explicit replacement behavior remain.

Use existing independent default-reference orientation candidates without a
synthetic study AF column. Extend preparation/row contracts to support absent AF
while preserving exact native-AF binding. Match the external source once per
chromosome shard with the existing complete-site allele join; retain complete
assignments for all possible final keys and all source rows, including rows
ineligible for diagnostics. Reuse the existing checked cache, atomic streaming
writer, source identities and worker adoption machinery rather than copy them.

## Scientific contract

Screen raw external counted-allele frequencies before any swap. Pool deduplicated
usable and AF-at-or-below-0.5 counts across observed study chromosomes, never
average fractions or count one shared source twice. Unobserved external
chromosomes do not belong to this dataset's screen. Zero usable frequencies fail.
If the pooled raw screen is not MAF-like, retain the existing declared ALT/EAF
contract and record that no independent semantic comparison was required.

For a MAF-like source, compare raw external counted-allele AF x against
independent default-reference AF r for that SAME external counted allele.
With assignments in default ALT direction and external assignment swap s:
`x = 1 - assigned_AF if s else assigned_AF` and
`r = 1 - default_ALT_AF if s else default_ALT_AF`.
The original study/default swap is not the external semantic direction.
Require unique, non-palindromic, biallelic eligibility in BOTH panels before
complementing or folding. Preserve complete reference-site membership before
filtering. Known multiallelic and unresolved sites cannot enter diagnostics,
but remain available to existing output attachment and rejection policies.

Reuse the existing dual-rule calculation and canonical policies:

| Requirement | Rule 1 | Rule 2 |
|---|---|---|
| Correlation | corr(x, r) >= 0.95 | corr(x, min(r, 1-r)) >= 0.95 |
| Mean absolute difference | mean(abs(x-r)) <= 0.10 | mean(abs(x-min(r,1-r))) <= 0.10 |
| Reference profile | proportion(r <= 0.5) > 0.90 | proportion(r > 0.5) > 0.20 |

The shared gate is at least 1,000 eligible finite pairs and pair count divided
by ALL study rows entering this preparation >= 0.50. Keep that fixed denominator
through exclusions. Run neither metrics nor ranking when this gate fails.
Evaluate both complete rules on the same pooled evidence. Rule 1 alone permits
the declared external source. Rule 2 alone stops with external MAF support;
there is no original study AF for the explicit replacement contract. Both pass
or both fail stops as inconclusive. Never auto-select the default or another
file, infer replacement, or relax thresholds to make the EGG data pass.

An independent reference is required for a MAF-like source: the same physical
file, including aliases or a different population column, is not independent.
The old 0.02 advantage margin and old chromosome-specific semantic decision do
not participate in approved dataset workers. Direct helpers without a dataset
approval retain their documented behavior. These thresholds are configured QC
heuristics, not proof of column semantics or empirical calibration.

## Scale, adoption and failure

Process narrow chromosome partitions one at a time. Spool replayable paired
evidence; use the existing stable pooled Pearson statistics, not averages of
chromosome correlations or an eager all-GWAS concatenation. Exact Spearman must
reuse existing global rank budget checks after count/coverage and may collect
only within that contract. No new hidden eager fallback or capacity claim.

Keep candidate caches private until global acceptance and unchanged-input checks.
Publish prepared default matches and complete assignments atomically into the
existing typed work shards. Validate the OLD source partition identity immediately
before each atomic replacement; update fingerprints afterward. Recheck all
bindings and approve only after every publication succeeds. Partial publication
never authorizes workers. Original/staged external files, declared aliases,
default files, mappings, policies, software, row identities, descriptors and
cached values remain bound. Missing/stale/corrupt approval must fail without
source rereading or falling back to a weaker classifier.

Workers validate approval before independently resolving strand, consume only
the assignment selected by the final key, and retain all current effect/strand,
palindrome, complete-site REF, duplicate, per-chromosome usable coverage,
missing/invalid/degenerate AF, rejection and output concordance checks. Do not
repeat semantic classification or MAF warnings for globally accepted AF.

Preflight failure raises the existing dataset-level `StudyAlleleFrequencyError`
with measured evidence, finalizes canonical failure reporting and starts zero
workers/retries. Extend current terminal, hidden-screen, manifest and HTML
reporting from recorded evidence only. No new threshold/configuration key is
needed; update affected canonical YAML descriptions and generated policy docs.
Do not edit protected adapters or unrelated user files/copies.

## Verification and known limits

Use RED/GREEN service-level missing-AF regressions, genuine-MAF and inconclusive
controls, independent source/default swaps, pooled count and coverage boundaries,
same-file aliases, native/replacement invariance, multiallelic/palindrome
eligibility, complete output/schema/QC parity, private-name collisions, subsets,
source/cache/config/partition mutations, partial publication, logging and
zero-worker/no-retry failures. Approved workers must do zero new external
parser/heavy matching or semantic classifier calls.

The available real five-chromosome EGG extract has N=2,077,782 and 972,784 usable
pairs (0.46818386144455965). It must remain inconclusive under the approved 0.50
gate. Its direct correlation 0.9972718912258661, direct MAE
0.008319938630348566 and same-allele minor profile 0.942543257290416 satisfy
Rule 1's numerical/profile conditions separately. It is not a full-GWAS positive
test. Verify the fixed failure stage on actual data and use representative
synthetic positive fixtures. Full 22-chromosome/server execution and 100–200M
capacity remain unverified unless actually run.

Primary contracts: [GWAS-SSF](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex),
[IGSR ALT frequencies](https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/),
[VCF AF cardinality](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
These support allele definitions, not calibration of the retained thresholds.


## Verified pre-commit outcome — 2026-11-05

The permanent [validation snapshot](../plans/2026-11-05-whole-study-external-af.md#pre-commit-validation-snapshot--2026-11-05)
contains exact commands, all 14 owned paths, hash convention, all seven rulings
and costs, independent row/numeric/QC oracles, scientific sources, tests and limits.
Final whole-task review, exact-path local commit and owned scratch cleanup remain
pending at that snapshot; later Git/final handoff supply those statuses.

Measured 403-source SHA256:
`8bb78f6f31e65f8264baf1c0811c590a9b7d25f104b7be3e71e58ba288aafdfa`;
source plus focused test:
`49552b3c1d6ad67acd0be9645c4c83d7ef57ac48c645b59a3a9e136fc55eb200`.
No configuration keys/defaults/constraints or protected adapters changed;
canonical YAML differs only in `help`/`origin` descriptions.

The actual five-chromosome final global entry correctly stopped at
`insufficient_overlap_fraction`, N=2,077,782, n=972,784, .46818386144455965,
with zero metric/legacy/AF-worker/strand-worker calls, no approval/replacement and
unchanged source/payload/shard bytes. The old margin reason was absent. The direct
entry harness and separate zero-worker/no-retry service regressions have distinct
scopes; this expected failure does not establish successful source acceptance.

Independent 4K/100K/1M two-chromosome checks passed exact frames/schema/order/dtypes,
row and pooled numeric oracles, explicit QC accounting and zero worker rereads.
Total external matches remain 2→2. At 1M, AF phases took 6.058029→18.228685s,
sampled RSS 987,430,912→2,931,163,136 bytes, with 4,638,717 added prepared-partition
bytes. These are functional preflight costs, not an optimization or capacity proof.
Both independently reproduced P2 findings were fixed and scoped re-review passed.
Final covering: 284 passed/131.79s/no warnings; broad 34 owners: 1,185 passed/
159 warnings/460.88s (unchanged warning producers); clean portable 10 owners:
492 passed/266.01s/no warnings. AST/YAML/policy/50-page wiki/protected-path checks
passed. Unrestricted retains exactly the 37 baseline collection-error paths,
exit 2; CLI extraction retains exactly two existing installer false positives,
exit 1. No full-suite release-green claim is supported.

Full22/all-11 workflows, full native production pipeline, current cloud/server
identity, server/Linux behavior and 100–200M capacity remain unverified. Exact
Spearman retains the existing rank-budget gate; sampled RSS is not a hard bound.
Primary contracts above support allele semantics, not threshold calibration.

Accuracy confidence: **98%** for recorded scoped behavior, bounded by the explicit
unperformed workflow/capacity checks and disclosed baseline tooling defects.

## Final wording correction gate — 2026-11-05

The two P3 findings in the independent whole-task review are corrected: the
shared replay-count diagnostic names `AF replay`, and the `_prepare_only`
docstring identifies both initial-external and explicit-replacement preparation.
No scientific decision, numeric path, default or configuration changed. The
measured source and tests above retain their exact `8bb78f6f...` attribution;
they were not rerun on the later wording bytes. The focused final gate passed
seven cases, and exact two-string AST normalization proved the runtime structure
unchanged. The final 403-source SHA256 is
`e06b266b7e06cccf14aebba576660ab12bba4a37d34935677d6feb6d5b54e6c5`;
source plus focused test is
`c7024ab5e25c5a8f00baa919850bd6797e29b3b25b6db3d32edc2bc16d78e1b5`.
The private `final-fix-source.json` gives component-wise `Path`-ordered file
identities; `final-fix-report.md` gives the exact focused command and result.
Final scoped review, local commit and scratch cleanup remain pending.
