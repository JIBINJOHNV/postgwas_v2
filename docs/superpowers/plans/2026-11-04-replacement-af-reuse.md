# Replacement AF Assignment Reuse Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Eliminate the second external replacement read/match while preserving
scientific outputs and validation decisions.

**Architecture:** Prepare complete row/key-bound assignments once. After the
global numeric gates pass, publish them atomically in existing chromosome work
partitions with checked input identities; issue approval only after all writes
and final identity checks succeed. Select the cached assignment after strand
processing and share existing finalization.

**Tech Stack:** Existing Python, Polars, Parquet, canonical YAML and pytest.

**Spec:** docs/superpowers/specs/2026-11-04-replacement-af-reuse-design.md

## Global Constraints

- No scientific cutoff, method, source-selection or public CLI/default change.
- Original AF and original diagnostic fixed N/eligible pairs remain unchanged.
- No replacement AF may resolve strand/palindromes; effects transform once.
- Preserve complete multiallelic REF sums, reasons and value-state handling.
- Recompute worker row/coverage statistics on its actual surviving subset.
- Cache all potential output keys, not just eligible diagnostic pairs.
- Reject stale/altered/incomplete caches without reread fallback.
- Reuse existing merge/statistics/atomic-write mechanisms; no duplicated logic.
- Preserve user files and protected adapters; no production workflow or push.

## Review Focus

1. Diagnostics exclude rows that workers may retain; complete cached output needed.
2. Upstream filters change worker denominators and unresolved-site reason counts.
3. Cache publication changes work fingerprints and can fail partway through.
4. Source aliases or in-memory cached values can change after measured approval.
5. Appending nested cache columns must not eagerly load a wide full study or
   change original row order/schema/content.

---

### Task 1: Complete assignment reuse, integration and regressions

**Files:** modify allele_frequency.py and shared/allele_join.py under
src/postgwas/modules/harmonisation; extend chromosome_partition.py only for
shared checked streaming atomic publication. Add no parallel implementation.
Extend tests/test_harmonisation_replacement_af_integration.py and
tests/test_harmonisation_frequency_reference_reuse.py; add focused atomic-writer
cases to the existing integration owner. Modify other owning files only when a
traced dependency requires it. Update affected module documentation/config help
only where existing behavior descriptions change; no default/schema change.

**Interfaces:** retain public merge behavior and return schema; generalize shared
read/join preparation and post-join finalization so prepared assignments can be
consumed without a second source read. A private row/key-bound cache descriptor
and its updated partition fingerprint belong in the existing replacement
approval binding. Workers consume that descriptor only under passed validation.

- [x] Read owning files/call paths completely; inspect AGENTS/state/diff.
- [x] Add failing end-to-end tests for one total distinct read/join and zero
  worker reads, complete eligible/ineligible assignment reuse, subset counts,
  integrity and atomic failure behavior. Verify expected RED on unchanged code.
- [x] Generalize existing preparation/finalization and raw joined-assignment
  metadata; share scientifically meaningful counts rather than copy them.
- [x] Derive diagnostics from the complete assignment table with the same n/N,
  original-EA B and metrics. Stage narrow caches privately; after global success
  and unchanged-input checks, append them atomically and update bindings.
- [x] Select cached assignments after existing strand processing. Require exact
  known row/key entries and run shared finalization on surviving rows. Keep
  identical-source reuse; remove all private columns from outputs.
- [x] Verify source/cache identity at approval and adoption, including alias
  mutation. No silent fallback and no worker approval after partial publication.
- [x] Run focused numerical/integration/reuse/multiallelic/downstream/reporting,
  partition and required-resource tests with the analysis bin on PATH.
- [x] Update affected docs and finish cumulative scientific/duplication review.
  Save exact commands/results and remaining limits for independent review.

### Task 2: Independent verification, benchmark and final handoff

**Files:** plan/spec verification record and affected module docs; root owns
verification artifacts. Runtime files stay frozen unless a reviewed defect needs
the original implementer to correct it.

**Interfaces:** consume Task 1's complete tested patch, cache protocol and report.
Produce baseline/optimized equality and source-I/O evidence, scientific review,
final scoped checks, limitations and a local follow-up commit when authorized.

- [x] Review the complete patch against Task 1 requirements independently.
- [x] Benchmark generated preparation plus adoption, including one distinct source
  with direct/swapped and multiallelic rows; compare full outputs/scientific QC
  with baseline. Record elapsed time, sampled RSS, file/cache bytes and reads.
- [x] Run broader tracked-plus-feature tests and unrestricted pytest; preserve
  baseline failures by exact name rather than repair unrelated code.
- [x] Verify syntax/configuration/docs, exact final diff and protected paths.
- [x] Complete independent whole-task review and address material findings.
- [ ] Finish root-owned exact-path local commit and task-created scratch cleanup;
  report the resulting status without implying a push.

## Permanent implementation and validation record (2026-11-05)

This is a **PRE-COMMIT VALIDATION SNAPSHOT**. Pending root actions below describe
the state at validation; subsequent commit and cleanup status is recorded by
Git history and the final handoff.

Implementation, whole-task review, scoped correction re-review and the final V4
benchmark are complete against integration base
`2e46f878eaa545515fe2aa4ac0ec6f97bf0e644b`. Final scoped regression/static
checks passed. Exact-path local commit and scratch cleanup remain; no follow-up
commit ID is invented. Repository-wide baseline blockers are disclosed below; no push or production execution is implied.

### Owned scope and preserved behavior

The exact ten owned paths are:

```text
docs/modules/harmonisation/configuration.md
docs/modules/harmonisation/eaf_replacement.md
src/postgwas/modules/harmonisation/allele_frequency.py
src/postgwas/modules/harmonisation/chromosome_partition.py
src/postgwas/modules/harmonisation/shared/allele_join.py
tests/test_harmonisation_eaf_replacement_science.py
tests/test_harmonisation_frequency_reference_reuse.py
tests/test_harmonisation_replacement_af_integration.py
docs/superpowers/plans/2026-11-04-replacement-af-reuse.md
docs/superpowers/specs/2026-11-04-replacement-af-reuse-design.md
```

Complete source-row/final-key assignments are prepared once per distinct
chromosome replacement source, covering diagnostic-ineligible output rows too.
After global gates pass, checked streaming atomic publication appends the private
cache to the existing work partitions and updates their fingerprints. Approval
is issued only after all publications and final identity checks succeed. Declared
aliases, physical identities, descriptors and exact row-bound values remain
checked. Workers select the key established by existing strand processing and
use shared finalization on surviving rows, with no replacement-source reread or
silent fallback. Identical-source reuse still requires zero additional reads.
The exact List<Struct> guard compares validity, boundaries, element order and
primitive fields losslessly; matching ordered source IDs skip remapping only.
Subsets/reordering retain their row-bound mapping and all integrity checks.

Independent review found Important I1: a mapped original-A column using a
private assignment name could be read as candidate B and incorrectly pass the
gate. Shared descriptor-driven name reservation now protects actual receiving
columns in preparation and adoption. The 22 collision cases cover valid and
constant-0.25 sources across eleven mapped names; the covering run passed 161
cases in 89.81 s. Independent re-review passed specification and code quality,
with no remaining scoped findings. Constant replacement sources issue no approval
or cache publication; original A, fixed N, metrics and cached/fresh outputs agree.

Final whole-task review then identified F1/F2/F3 and prose correction M1. F1 now
revalidates the old partition/source binding after streaming/output checks and
immediately before atomic replacement, so a scientific input mutation is refused
without overwriting the independently modified file. This is not a general
cross-process lock or crash-durability guarantee. F2 allocates the public AF
destination against scientific input before reserving private assignments and
markers, preserving literal fresh-output schemas/mappings for configured names
and suffix chains. F3 passes actual matching orientations through shared
finalization and removes preparation-only site/reconstruction QC for direct-only
attachment; diagnostic eligibility remains complete. M1 corrects the Architecture
approval order above. All four findings are addressed; independent scoped
re-review passes specification and quality with no remaining findings. The
23 new regressions pass in 25.06 s and 253 covering cases pass in 117.48 s,
without a warning summary. Existing I1 cases remain covered; assertions and the
two permitted I/O differences were not weakened.

No scientific cutoff, policy, source-selection rule, configuration key/default,
CLI or dependency changed. Protected adapters remain untouched. The six modified
Python owners compile in memory. AST comparison preserves
`validate_study_allele_frequency`, `_measure_frequency_pairs`,
`validate_replacement_frequency_pairs`, `_compare_frequency_pairs`,
`_reference_frequency_sites` and `harmonise_strand_orientation` unchanged.
Scientific sources and interpretation limits remain in the
[replacement guide](../../modules/harmonisation/eaf_replacement.md); streaming
behavior is supported by the versioned official Polars reference linked there
and direct installed-0.20.31 tests, rather than a new scientific method.

### Authoritative V4 benchmark

Baseline Git source and the final frozen 403-file source copy used the same
physical generated Parquet inputs under `/private/tmp`, each in a fresh process.
Python 3.10.21, Polars 0.20.31 and NumPy 1.23.5 ran on Darwin arm64, 12 logical
CPUs/64 GiB RAM, with four Polars threads. RSS was sampled every 20 ms.
Final `allele_frequency.py` SHA-256:
`f356f0030c50bc0c72f6bea1e7894576b8783bceb03332187d6aec6fbd9cc723`.

| Study rows | Baseline s | Final V4 s | Less elapsed time | Baseline sampled RSS B | Final sampled RSS B | Baseline partition B | Final partition B |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100000 | 3.5297437501 | 2.9340127502 | 16.88% | 717570048 | 629686272 | 890661 | 1598671 |
| 1000000 | 32.5472067501 | 20.0146634171 | 38.51% | 4743495680 | 3930423296 | 8415982 | 14543715 |

Extra compressed partition storage is 708,010/6,127,733 bytes; this is not peak
temporary disk use. Sampled RSS is 12.2%/17.1% lower, not measured OS peak RSS.
The workload uses 110K/1.1M replacement rows, all four direct/swapped/complement
orientations and 10% multiallelic sites requiring complete-site REF reconstruction.
Timing includes original AF validation, initial partition writing, replacement
preparation/publication and adoption; it excludes input generation, initial study
loading, imports and proof serialization. All rows survive and 90% enter diagnostics;
ordered-ID performance is measured, while subset/reordering performance is not.

Generation is deterministic: for `i=0..N-1`, chromosome is `1`, source ID is `i`,
`POS=100+(37*i)%N`, `MAF=.05+.4*(i%1000)/999`, and default ALT AF is
`1-MAF` when `i%3==0`, otherwise `MAF`. Original EA/OA cycles through
`C/A,A/C,G/T,T/G`; `BETA=.1+.001*(i%23)` and `Z=BETA/.05`. Default REF/ALT is
`A/C`. For `i%10==0`, replacement REF/ALT is `C/A` with AF
`1-default_AF-.05`, plus an extra `C/G` row at AF `.05`; other rows equal default.
Inputs use Zstd Parquet; references sort by POS and replacement additionally by
ALT. Canonical `default_policies()` are unchanged. Column mapping is
CHROM/POS/ALT/REF/AF; original MAF, beta effect and forward-strand metadata feed
`validate_study_allele_frequency`, a Zstd work-partition write,
`validate_replacement_allele_frequency`, then `harmonise_allele_frequency` on the
reloaded partition. Scratch benchmark-script invocations are historical task
commands, not installed/public reproducers; this recipe, code and tests preserve
the experiment's reconstruction basis after task-workspace cleanup.

All complete 13-column output frames match old code exactly in schema, order and
values, with no numerical tolerance or observed roundoff differences. Scientific
QC agrees exactly except the declared I/O fields:
`external_merge_external_rows_loaded` becomes zero in workers and
`external_merge_reused_prepared_replacement_assignments` becomes true. Comparison
excludes the implementation/workfile-bound approval contract and normalizes only
the per-run duplicate-report basename; other identity/scientific fields agree.
Replacement parses, merges and heavy joins each fall from two total to one,
with zero worker additions. Three independent old-code differentials also pass
exactly: missing multiallelic component, invalid ALT sum and retained
default-unmatched row; each retains 12 rows and 11 diagnostic pairs.

The first candidate was rejected: although equality and read counts passed,
time increased about 60% and sampled RSS 2.35–2.77-fold. Profiling identified
native nested equality in the approval guard; the lossless comparison corrected
that bottleneck. V3 preceded the final publication-identity, destination-name
and direct-only-QC corrections; V1/V2/V3 are historical, not final performance
evidence. V4 measures the fully corrected frozen source. An initial Documents-
folder 1M baseline stopped on identity drift of unknown cause; controlled inputs remained
stable and checks were not weakened. These single-run local AF-path observations
do not establish full harmonisation, production/Linux or 100–200M capacity,
cutoff calibration, general speed/RAM guarantees or per-row EAF semantic proof.

### Regression, static checks and remaining limitations

| Check | Result |
|---|---|
| Isolated baseline, nine affected owners from Git HEAD | 317 passed in 106.30 s |
| Final broad workspace, 32 owners | 1079 passed, 159 warnings in 413.76 s; exit 0 |
| Exact clean candidate, HEAD plus ten owned paths, nine affected owners | 421 passed in 203.36 s; exit 0; 104 new cases over baseline |
| Unrestricted pytest and exact baseline error-set comparison | 37 existing collection errors in 12.49 s; exit 2; exact paths unchanged; one third-party warning |
| Final syntax/AST/configuration/protected-path/doc checks | Passed; six Python owners, six unchanged scientific ASTs, policy synchronization and 50-page wiki build |
| Independent whole-task/correction reviews | F1/F2/F3/M1 addressed; specification and quality PASS; no remaining findings |
| Local commit and scratch cleanup | Root-owned final repository operations; pending |

Root invocations use analysis `bin` first on PATH, candidate-only `src` for the
isolated run, `PYTHONDONTWRITEBYTECODE=1`, Python `-B`, and
`pytest -q -p no:cacheprovider` with owned temporary `--basetemp` directories.
All modified runtime/test bytes match the tested clean candidate; all 403
benchmarked source files still match current source. The 159 warnings originate
at nine unchanged producer locations: seven concordance left-join deprecations
(`analysis.py:205,814,822,1071,1096,1624,1641`) and filename/file-object warnings
(`summary_statistics_io.py:208`, `test_harmonisation_metadata_preamble.py:178`).
No warning producer is in touched code; the count is not claimed equal to a
differently scoped baseline.

Equivalent completed commands (the actual runs additionally used owned
`--basetemp` paths), with the analysis environment first on PATH and bytecode
disabled as above. The candidate contains only Git HEAD plus the ten owned paths:

```bash
# broad: run from repository root
PYTHONPATH=src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_delimiter_header_selection.py tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_build_fallback.py tests/test_harmonisation_config.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_exclude_chr_x.py tests/test_harmonisation_exclude_x_accounting.py tests/test_harmonisation_exclude_x_interface.py tests/test_harmonisation_exclude_x_references.py tests/test_harmonisation_exponentiated_se.py tests/test_harmonisation_exponentiated_se_integration.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_metadata_preamble.py tests/test_harmonisation_multiallelic_reference_frequency.py tests/test_harmonisation_reference_frequency_export.py tests/test_harmonisation_reference_frequency_provenance.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_resource_preflight.py tests/test_harmonisation_sample_sheet_generator.py tests/test_harmonisation_shared_coordinate_identifier.py tests/test_harmonisation_study_report_integration.py tests/test_harmonisation_allele_join.py tests/test_harmonisation_imputation_quality.py tests/test_harmonisation_eaf_replacement_downstream.py tests/test_harmonisation_chromosome_partition_index.py tests/test_harmonisation_strand.py tests/test_harmonisation_external_reference_staging.py tests/test_harmonisation_adaptive_scheduling.py
# candidate: run from the isolated candidate root
PYTHONPATH=src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_resource_preflight.py
# unrestricted: run from repository root
PYTHONPATH=src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider
```

Evidence artifacts are `root-final-verification-commands.json`, its results/logs,
`task-1-report.md`, `task-1-review.md`, `task-1-rereview.md`,
`final-fix-report.md`, `final-fix-rereview.md`,
`final-optimized-v4-comparison-summary.json`, `edge-case-v4-comparison.json`,
`optimized-source-checksums.json` and `final-benchmark-record.md` in this plan's
task workspace. Their material outcomes are preserved here before cleanup.

Unrestricted collection retains all 37 baseline error paths, none added or
removed. The one unrestricted warning is the third-party marshmallow distutils
deprecation. Exact error paths and the installer-example extractor defect are in
the [preceding validation plan](2026-11-04-replacement-af-validation.md#explicit-unresolved-baseline-checks).
The fresh CLI checker again exits 1 on the same two installer false positives in
`README.md` and `docs/wiki/getting-started/installation.md`: its extractor reads
`postgwas --jobs 2` from valid installer `--name postgwas --jobs 2` arguments.
Root rechecked the extractor and installer support for `--jobs`; neither file
changed. The local broad run also uses an existing ignored concordance helper
and a partition-index owner referencing private frozen history; it is not a
whole-suite fresh-clone proof.
The isolated nine-owner candidate specifically tests affected portable owners.
No unrelated repair, fixture payload or ignore-policy change is included.

### All five ledger rulings and costs

1. Ruling: Reuse the active user feature checkout, preserving unrelated files — native workspace instructions prefer the current suitable checkout — wrong attribution requires correcting staged hunks.
2. Ruling: User's explicit optimization request approves the described complete-assignment reuse scope; routine cache/mechanism choices proceed without another approval round — no new scientific policy — a mistaken engineering choice costs visible rework, not a silent scope change.
3. Ruling: Embed checked cache in existing chromosome work shards using shared atomic writer and verified streaming append — avoids unmanaged persistent sidecar paths and new defaults — cache adds disk/streaming-join state, so memory and timings must be measured and bounded claims withheld.
4. Ruling: Existing unrestricted collection/doc-validator defects remain disclosed baseline limitations; in-scope regressions are fixed without an approval pause — unrelated fixes excluded — no full-suite release-success claim.
5. Ruling: Replace slow nested protocol equality with lossless flat comparison of schema, list validity/length, element validity and primitive fields; retain exact comparisons, no hashing/tolerances, plus safe same-ID-order fastpath — measured bottleneck fix rather than weakening checks or redesigning storage — if wrong invalid rows could be accepted/rejected, requiring native-oracle/null/boundary/corruption tests and full-output differential benchmark.

Compliance: scoped implementation/documentation checklist passed; cumulative
diff reviewed; hardcoding/duplication audit passed; scientific validation
completed for the recorded scope. Whole-task review and scoped correction
re-review are complete; root-owned local commit and scratch cleanup remain at
this pre-commit snapshot, with subsequent status in Git history/final handoff.
Unrestricted collection and the CLI checker remain blocked by the explicit baseline issues above. No full-suite success is
claimed. Documentation whitespace and local-link checks pass; no runtime tests
were repeated for prose-only bookkeeping. No new ruling or deferred minor is
introduced by this record.
