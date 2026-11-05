# Reuse complete replacement allele-frequency assignments

The user authorized optimizing the repeated work identified after commit
2e46f87. Retain complete replacement assignments during whole-study preparation
and consume them in chromosome workers. Scientific decisions, thresholds,
frequency values, effect transformations and row-retention policies must stay
unchanged. This extends the previously approved replacement design; no new
user preference or scientific rule is needed.

## Approach and scope

Cache complete, narrow row-bound assignments inside the existing typed chromosome
work partitions. This avoids a new persistent sidecar/output-layout policy and
lets existing prepared-work fingerprints, scheduling and checkpoints cover the
cache. Keeping only diagnostic pairs is insufficient: those pairs exclude
palindromes, multiallelic sites and unmatched rows that may survive workers.
Caching only a normalized reference would still repeat the expensive join.

Use all prepared default-reference candidates plus canonical forward and reverse
complement possibilities for default-unmatched rows. Deduplicate possible keys
by source-row ID and canonical chromosome/position/EA/OA. Replacement frequencies
must never select an orientation. Store final-key frequency, match orientation,
complete-site reason, biallelic eligibility and any multiallelic flag needed for
the existing counters; preserve static source/duplicate metadata separately.
Retain original-EA restoration for the unchanged diagnostic subset and fixed N.

Generalize the current merge into shared preparation and finalization rather
than duplicate its matching, frequency-site arithmetic, statistics or coverage
logic. Ordinary external sources retain the same public behavior. Workers use
the actual post-strand key and stable source ID to select exactly one cached
entry, including explicit unmatched/null assignments, then apply the shared
finalization on their surviving rows. This recomputes meaningful counts and
coverage after upstream exclusions. Remove private cache fields before outputs.

## Publication, integrity and resources

Keep candidate caches private until all existing global gates pass and original
bindings remain unchanged. Append caches atomically to existing work shards,
update their fingerprints, recheck sources/defaults/staging and issue approval
only after every publication succeeds. Include cached values in the existing
row-bound worker comparison. Missing, malformed, truncated, altered or stale
cache data must fail; never silently reread the replacement source.

Preserve source-alias checks, existing original-input evidence, legitimate row
subsets/reordering, multiallelic REF=1-sum(ALT), palindrome orientation requirements,
original MAF/BETA/Z provenance, and worker coverage failure behavior. A partial
cache publication cannot authorize workers. Reuse the existing atomic Parquet
writer; extend it for a checked streaming append when supported by the installed
Polars runtime, without an eager wide-study collect or hidden eager fallback.
Keep work one chromosome at a time; global Spearman resource limits are unchanged.
Cache storage and streaming join state remain real costs to measure.

## Verification

The observed baseline distinct-source path uses two external merge/parser calls
across preparation and adoption; the identical-source path uses zero. Require one
total distinct-source match/read per chromosome and zero in workers. Compare
actual frequency values, scientific QC counts, rejected rows, effects and output
order with baseline. Only explicit cache/I/O provenance may differ.

Cover direct, swapped and reverse-complement alleles; palindromes; retained
reference-unmatched rows; multiallelic direct and reconstructed REF; invalid or
incomplete sites; missing/invalid frequencies; source-column collisions;
subsets/reordering; source/partition/cache mutation; alias retargeting; failed
global gates/publication; and separate chromosome coverage. Benchmark the full
preparation-plus-adoption layer on generated data and record time, sampled RSS,
cache bytes and source-read counts. No full-workflow 100-200M capacity claim.

The [Polars sink_parquet documentation](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.LazyFrame.sink_parquet.html)
describes streaming writes. Installed Polars 0.20.31 was separately characterized
with a nested list/struct left join and streaming sink; schema and original row
order were preserved. Tests must exercise this installed-runtime path directly;
documentation alone does not establish compatibility or a RAM guarantee.


## Exact prepared-row comparison

Profiling the first candidate showed that native equality on the two nested
List-of-Struct protocol columns dominated worker time and allocation. Compare
these losslessly by schema, outer-list validity, per-row lengths, exploded
element validity and flat field equality. Preserve row and element order and
NaN/null semantics; never substitute a hash or numeric tolerance. When source IDs
match in physical order, skip unnecessary remapping/sorts; genuine subsets or
reordering keep the existing row-bound mapping. Native equality on small
structural/mutation cases and full old/new outputs are required oracles.

## Implemented verification boundary (2026-11-05)

Complete assignments now share the existing merge/finalization path and checked
streaming work-partition publication. The mapped original-A/private-name
collision identified in independent review is fixed by shared descriptor-based
name reservation; 22 collision cases and 161 covering cases pass, followed by a
clean independent correction review. Native-equality oracles cover the exact
List<Struct> guard; matching ordered IDs avoid remapping without omitting values.
Source/partition aliases, descriptor integrity and row-bound cache checks remain
mandatory. Final F1/F2/F3/M1 corrections protect the pre-replacement input
identity, literal public destination names, mode-specific direct-only QC and
approval-order wording. All four pass independent specification/quality re-review;
23 focused and 253 covering cases pass. No scientific policies, cutoffs, defaults,
CLI or dependencies change.

The authoritative V4 comparison has exact complete 13-column old/new frames and
scientific QC apart from two explicit worker I/O fields, one rather than two
replacement reads/joins, and exact independent special-case differentials. Local
100K/1M AF-path time falls 16.88%/38.51%, sampled RSS falls 12.2%/17.1%, and
compressed partitions add 708,010/6,127,733 bytes. These single-run observations
do not establish production/Linux or full-workflow 100–200M capacity. The
[permanent plan record](../plans/2026-11-04-replacement-af-reuse.md#permanent-implementation-and-validation-record-2026-11-05)
preserves exact measurements, all five ledger rulings, review outcomes, owned
paths, 1079-pass broad and 421-pass isolated-candidate results, and the unchanged
37-error unrestricted baseline. Final static, policy and wiki checks pass; the
known CLI-checker false positives remain disclosed. Whole-task review and scoped
correction re-review are complete. In this **PRE-COMMIT VALIDATION SNAPSHOT**,
local commit and scratch cleanup remain; subsequent status is recorded by Git
history and the final handoff.
