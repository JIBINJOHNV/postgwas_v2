# Whole-study Initial External AF Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Correct the no-mapped-study-AF path so its initially selected external
source is validated once across the whole dataset and reused by workers.

**Architecture:** Generalize existing optional-AF reference preparation and the
checked replacement assignment/cache machinery. Reuse pooled dual-rule metrics
for same-external-counted-allele evidence; publish only after global acceptance.
Keep native AF, explicit replacement and direct-helper compatibility unchanged.

**Tech Stack:** Existing Python, Polars, Parquet, canonical YAML and pytest.

**Spec:** docs/superpowers/specs/2026-11-05-whole-study-external-af-design.md

## Global Constraints

- No new thresholds/resources, hidden fallback or relaxed coverage; approved
  defaults remain 1000, 0.50, 0.95, 0.10, 0.90 and 0.20 in canonical YAML.
- Raw-source MAF screen precedes any swap; pool counts over observed chromosomes.
- Compare raw external ALT AF to independent reference AF for that SAME allele.
- Require both panels' complete-site biallelic eligibility before diagnostics.
- Keep fixed N, all scientific input rows, complete output candidates and effects.
- Rule1 alone accepts; rule2 alone stops; both/neither stops; no automatic source
  replacement when no native study AF exists.
- Reuse matching, metrics, cache, atomic writer, identities and reporting owners.
- Stream narrow chromosome work and Pearson batches; reuse Spearman resource gate.
- Never approve partial publication; recheck OLD partition immediately pre-replace.
- Preserve every native/replacement guard and worker QC; no legacy classifier
  fallback or second source read/match in approved normal dataset workers.
- Preserve user untracked copies/directories; protected adapters untouched.

## Review Focus

1. Initial external AF and explicit replacement are different modes; no fake
   original AF/evidence or borrowed replacement flag may authorize adoption.
2. External/default/original study swaps are independent; compare the same counted
   allele and protect both panels' complete-site eligibility.
3. Cache publication/mutation/partial failure and descriptor naming must preserve
   exact scientific output, resolved columns, provenance and source identities.
4. Multiple chromosome shards and shared staged files must not dilute/double-count
   raw screen or pooled metrics; no eager all-GWAS Pearson concat.
5. Deterministic global rejection must preserve evidence in logs/HTML/manifests
   and start zero workers/retries; the real extract's 46.8% gate must not be hidden.

---

### Task 1: Generalized initial-external preparation, adoption and regression tests

**Files:** modify existing owning functions in
`src/postgwas/modules/harmonisation/allele_frequency.py`, `strand.py`, `service.py`,
`study_reporting.py` and `html_report.py` as traced. Change `chromosome_partition.py`
or shared statistics/join owners only if necessary for focused reuse. Update
canonical `src/postgwas/config/defaults/modules/harmonisation.yaml` descriptions,
owning module docs and generated policy documentation with the established tools;
no new parameter keys/defaults. Add a portable focused test owner
`tests/test_harmonisation_external_af_preflight.py` and extend existing AF,
replacement, resource/reporting/service owners only where their interface changes.

**Interfaces:** produce a public preparation entry point for selected external AF
after staging; existing `validate_replacement_allele_frequency` remains a supported
replacement entry point/no-op for unrelated modes, delegating shared mechanisms
where appropriate. New separate `eaf_external_validation` evidence binds ordinary
approval. Use `resource_preflight['resources_by_chromosome']` to retain original
external paths after staging; do not add independent source resolution.
Generalize reference preparation to nullable native AF without synthetic AF data.
Generalize shared comparison to replayable batches with known count/fixed N and
resource-checked global ranks; existing DataFrame numeric behavior stays intact.

- [x] Read applicable AGENTS, actual owning files/call paths/config/tests in full;
  inspect state and cumulative diff; search existing mechanisms before helpers.
- [x] Write and run meaningful RED tests for absent mapped AF entering the global
  stage, pooled acceptance despite individual old-margin/count failures, exact
  raw external/default swap direction and no repeated worker parser/classifier.
  At least two chromosomes must jointly clear the count floor; fixed-N coverage
  equality passes, below 0.50 stops before metrics/worker calls.
- [x] Add FT/TT/FF and constant/missing/invalid controls, same physical source via
  alias/population column, complete-site/palindrome exclusions and unmatched rows.
  Keep zero-worker/no-retry assertions on real service error boundaries.
- [x] Generalize existing preparation/cache and streamed numeric helpers; keep
  normal native/replacement outputs and direct helper behavior unchanged.
- [x] Add approval-aware strand and AF adoption, same-source handling, checked
  source/partition/cache/descriptor/policy guards; retain subset behavior,
  collision-safe literal output schemas and all final per-variant QC.
- [x] Cover input mutation before/during publication, partial write failure and
  absent/changed approval with no source-reread/legacy fallback.
- [x] Integrate recorded evidence into service, terminal/hidden logs, manifests,
  HTML/failure summaries; update canonical descriptions/module docs/policy docs.
- [x] Run focused tests with representative fixtures; confirm RED/GREEN outputs,
  native/replacement invariance, policy/help synchronization and cumulative diff.
  Save exact commands/results, file scope and remaining limits; freeze runtime.

### Task 2: Independent validation and permanent evidence record

**Files:** plan/spec verification record and existing tests only if review finds
a material gap. Runtime frozen except reviewed corrections by Task1 implementer.

**Interfaces:** consume tested Task1 code/report; produce independently reviewed
current-code real-data and generated-scale output/evidence verification, broader
portable regression results and a concrete local commit for the authorized fix.

- [x] Independently review the full Task1 diff for spec/science/code quality;
  correct all material findings through the original implementer and scoped review.
- [x] Run current corrected global stage on the actual five-chromosome extract;
  require canonical insufficient-coverage failure with measured N/n, no workers,
  unchanged files and no old margin error. Do not claim full22/server success.
- [x] Run a representative generated multi-chromosome preparation/adoption
  comparison with fresh assignments: exact frames/QC, zero worker rereads and
  pooled evidence independent of chunk/chromosome order; record scale/limits.
- [x] Run affected portable owners and broader relevant owners, unrestricted
  pytest with exact baseline error accounting, syntax/config/doc checks; no
  unrelated collection/CLI-extractor repair or ignored-fixture promotion.
- [x] Preserve material commands/results/scientific decisions, unchanged defaults,
  all rulings/costs and explicit limitations in permanent plan/spec docs.
- [ ] Complete independent final review, exact-path commit and owned scratch
  cleanup; inspect final cumulative diff and preserve user files.


## PRE-COMMIT VALIDATION SNAPSHOT — 2026-11-05

This is the permanent evidence record before private task-workspace cleanup.
Task 1 implementation and its scoped correction review are complete. Task 2's
real-data, generated-scale, regression and static checks below are recorded.
**Final whole-task review, exact-path local commit and owned scratch cleanup are
still pending at this snapshot.** Subsequent final handoff and Git history supply
those later statuses; this document does not claim they already happened.
The private reports/scripts named below are provenance for executed checks and
are not durable dependencies or promised post-cleanup runnable files.

### Tested identity and exact file scope

Base: `1ea6305d54ca1343b9a871115018b0fa080396d2`;
branch: `codex/whole-study-external-af`.
Frozen measured scientific source: **403 regular `src` files**, SHA256
`8bb78f6f31e65f8264baf1c0811c590a9b7d25f104b7be3e71e58ba288aafdfa`.
The same calculation over `src` plus the focused test gives
`49552b3c1d6ad67acd0be9645c4c83d7ef57ac48c645b59a3a9e136fc55eb200`.
Convention: sort repository-relative `Path` objects by path components, excluding `__pycache__`; feed each
path's UTF-8 bytes, NUL, exact file bytes, NUL into SHA256. These are source/test
content identities, not commit identities or hashes of this self-changing record.
Final measurements use the corrected `8bb78f6f...` freeze; earlier `d6a2df36...`
and `62aa5368...` measurements are historical and do not support final timing claims.

All 14 owned paths (only the new focused test is force-add eligible under the
existing `/tests/` ignore rule; unrelated ignored fixtures are not promoted):

- `docs/modules/harmonisation/configuration.md`
- `docs/modules/harmonisation/eaf_replacement.md`
- `docs/wiki/harmonisation/policies.md`
- `src/postgwas/config/defaults/modules/harmonisation.yaml`
- `src/postgwas/modules/harmonisation/allele_frequency.py`
- `src/postgwas/modules/harmonisation/external_reference_staging.py`
- `src/postgwas/modules/harmonisation/html_report.py`
- `src/postgwas/modules/harmonisation/service.py`
- `src/postgwas/modules/harmonisation/shared/allele_join.py`
- `src/postgwas/modules/harmonisation/strand.py`
- `src/postgwas/modules/harmonisation/study_reporting.py`
- `tests/test_harmonisation_external_af_preflight.py`
- `docs/superpowers/plans/2026-11-05-whole-study-external-af.md`
- `docs/superpowers/specs/2026-11-05-whole-study-external-af-design.md`

The implementation reuses the existing reference preparation, complete-site
join/finalization, candidate assignments, cache/descriptor integrity, atomic
writer and input identity machinery. New public
`validate_external_allele_frequency()` runs after staging and before workers;
`capture_external_frequency_sources()` binds original aliases before staging.
Existing `validate_replacement_allele_frequency()` remains the separate supported
replacement entry/no-op in other modes. No native AF is synthesized. The distinct
`eaf_external_source`/`eaf_external_validation` metadata records reference-proxy
provenance and ordinary approval. Workers require current bound approval before
strand/AF adoption; incomplete or changed approval cannot trigger a weaker classifier.

Configuration keys added/removed: **none**. Recursive comparison against base
with only `help` and `origin` omitted proves all remaining YAML content identical,
including defaults and constraints. The scientific requirements remain 1,000
pairs, fixed-N fraction >=0.50, r>=0.95, MAE<=0.10, minor fraction >0.90 and major
fraction >0.20. The raw MAF-like screen and supported direct-helper controls are
also unchanged. No protected adapter, installer, unrelated test, public flag,
resource selection, automatic replacement or execution-default change is included.

### Scientific decisions and independent review

Raw deduplicated source counts are pooled over observed study chromosomes before
any swap. Diagnostic x and r both count the external source's declared allele:
external assignment swap complements both; study/default swap is independent.
Both panels must satisfy complete-site biallelic, unique and nonpalindromic
eligibility before arithmetic. All original study rows remain in fixed N and
complete output candidates remain available to normal final variant QC.
Rule 1 alone accepts; Rule 2 alone, both rules or neither rule stops. A MAF-like
source requires a physically distinct default file, including aliases/population
columns. A non-MAF-like source retains its declared ALT contract without an
independent semantic comparison. Thresholds are configured QC heuristics, not
empirically calibrated proof of frequency semantics.

Primary contracts were checked in Task 1 and independent review:
[GWAS-SSF v1.0.2](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex)
(effect-allele frequency),
[IGSR ALT frequencies](https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/)
and [VCF 4.5 AF cardinality](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
The installed GitHub connector independently retrieved GWAS-SSF blob
`c2b93ddc9804f86668dcf67e4af9fc515be5d32b`. These sources support allele definitions
and multiallelic handling, **not calibration of these thresholds**. Actual streaming
publication is exercised using Polars 0.20.31; its
[versioned sink contract](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.LazyFrame.sink_parquet.html)
was also checked.

The independent full Task 1 review found two P2 issues, both reproduced before
fixing and closed by scoped re-review: (F1) resource/typed failures lost measured
external evidence, and (F2) the previous chromosome's packed proof survived into
next preparation. The final fix preserves known N/source/count/coverage/resource
facts, null unmeasured metrics and unevaluated rules, including explicitly partial
raw counts on second-shard failure; it releases packed proof after its last spool
use on both source routes. No scientific threshold/calculation changed. Review
covered the entire original 14-file package and the complete three-file correction.
No material finding remains open in that scoped review; whole-task review is pending.

### Focused RED/GREEN and normal/boundary/failure coverage

Runtime for pytest: `/Users/JJOHN41/miniconda3/envs/postgwas/bin/python`, Python
3.10.21, Polars 0.20.31. Analysis `bin` is first on PATH; bytecode/cache writing is
disabled. No installation, version switch or canonical threshold override occurred.
Task 1's initial RED was 7 failures, followed by 7 GREEN; expanded controls exposed
and corrected native/replacement metadata, schema, service-entry, source-binding
and staging-alias regressions. Historical 332-test covering (192.90s) and 98-test
staging/alias check (65.62s) passed without warnings before the final formal fixes;
they are not substituted for the final frozen covering below.

Final implementer covering, inspected and confirmed by root:
**284 passed in 131.79s, no warnings**, exit 0:

```bash
PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src POLARS_MAX_THREADS=4 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q --tb=short -p no:cacheprovider --basetemp=/private/tmp/external-af-fix1-final-cover tests/test_harmonisation_external_af_preflight.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_af_acceptance.py
```

Formal-finding reproduction used the same environment/prefix, substituting
`--basetemp=/private/tmp/external-af-fix1-red` and these exact nodes:

```text
tests/test_harmonisation_external_af_preflight.py::test_spearman_resource_failure_preserves_measured_external_evidence_and_reports
tests/test_harmonisation_external_af_preflight.py::test_typed_initial_preparation_failure_finalizes_known_evidence
tests/test_harmonisation_external_af_preflight.py::test_previous_chromosome_packed_proof_is_released_before_next_preparation
tests/test_harmonisation_external_af_preflight.py::test_real_service_missing_native_af_enters_global_preparation_and_fails_before_workers[spearman_budget]
```

Observed six RED failures in 3.35s, then six GREEN in 4.38s with
`--basetemp=/private/tmp/external-af-fix1-green2`. The additional exact node
`tests/test_harmonisation_external_af_preflight.py::test_typed_second_shard_failure_retains_partial_raw_counts_without_claiming_global_screen`
failed in 1.27s with `--basetemp=/private/tmp/external-af-fix1-partial-red` before
correction, then joined the final 284-test covering. The .001-GiB Spearman case
retains N=n=1,200, coverage=1, raw counts=1,200, source identities, shared gate=True,
rules evaluated=False, null numeric metrics, estimated/available resource facts,
invalid status and no contract. The second-shard case retains only chromosome 1's
600/600 raw counts and labels preparation partial.

The focused owner and existing compatibility owners exercise:

- Two 500-row chromosomes jointly meet the 1,000 count equality boundary;
  1,200/2,400 coverage equality passes and 1,200/2,402 fails before metrics.
  Original fixed N cannot be reduced by partition metadata or exclusions.
- Independent external/study swaps, pooled acceptance despite the old .02-margin
  failure, Pearson replayed batches versus DataFrame and supported exact global
  Spearman including resource rejection. FT/FF controls stop; TT uses explicit
  profile/screen test overrides to make both profiles satisfiable and must stop.
- Same-file aliases and population columns; zero usable, missing, invalid,
  constant/degenerate and endpoint frequencies; shared-source observed chromosome
  scope/dedup; non-MAF-like exact-source reuse without extra external parse/join.
- Thirteen exact cached-versus-fresh frame/schema/column/QC comparisons: ordinary,
  direct-only, subset/reorder, private source name, suffix collision, missing,
  invalid, degenerate, palindrome, complete/unresolved multiallelic, unmatched and
  default-multiallelic cases. Only enumerated cache I/O accounting may differ.
- Missing approval and changed source/partition/mapping/column/descriptor/cache/
  approval/policy; mutation after metrics, after sink before replace and partial
  publication. The OLD partition is rechecked immediately before atomic replace;
  all publication must succeed before approval.
- Real service/resource-preflight/staging controls for passed/inconclusive/coverage,
  same-source success/failure, original content/alias, default alias, staged values/
  staged alias and mutation during staging; deterministic failures issue invalid
  evidence, failed stage status, zero workers and zero retry entry calls.
- Canonical terminal, visible/hidden screen logs, manifests and HTML use recorded
  external evidence without reanalysis. Native-AF/replacement compatibility owners
  remain unchanged; direct library helper behavior remains supported.

### Actual five-chromosome global-entry safety regression

Input: `/Users/JJOHN41/Downloads/gwas_sumstat/harmonisation_log_review_20261104/af_file_comparison/wlm_fetal_shared_frequency_fields.tsv.gz`,
31,435,578 bytes, SHA256
`b506002f53b623ef49a38c8817e5d06b76adeb4ed4a76d42542cdbc664cde85c`.
All eight columns were retained; `eaf_col=None`, so native `AF_INPUT` fields
remain unmapped. External 1000G EUR and independent panUKB EUR use the ten local
`GRCh37_{1000G,panukb}_freq_chr{6,7,19,20,22}.tsv.gz` resources in
`/Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/GRCh37/default_af/tab_files/`.
All eleven local payload fingerprints were verified before and after execution.

| Chromosome | Fixed study rows | Eligible pairs |
|---|---:|---:|
| 6 | 747,568 | 357,021 |
| 7 | 667,825 | 309,993 |
| 19 | 250,737 | 113,481 |
| 20 | 257,290 | 121,369 |
| 22 | 154,362 | 70,920 |
| Pooled | 2,077,782 | 972,784 |

The baseline actual chr22 helper at base commit reproduced
`error_difference_below_required_margin` with 70,920 pairs. Its folded r was
0.9953734403623264, folded MAE 0.008097881823011844, direct MAE
0.008190642398307954 and external-to-reference-MAF MAE 0.02646978565831923;
the canonical baseline margin remained .02 (current correlation minimum .90;
historical failure log .80). Separate baseline pooled diagnostic measurements gave
direct r=.9972718912258661, MAE=.008319938630348566 and same-allele minor profile
.942543257290416. These conditional numerical criteria do not establish acceptance.

The final corrected actual global entry stopped exactly as required:
`StudyAlleleFrequencyError`, `inconclusive/insufficient_overlap_fraction`,
972,784/2,077,782 = **0.46818386144455965**, below unchanged .50. At fixed N the
required count is 1,038,891, a gap of 66,107 pairs. Recorded exclusions reconcile:
449,858 original-orientation-ineligible + 655,140 external non-biallelic/unresolved
+ 972,784 comparable = 2,077,782. Raw reference screening has 19,736,516 usable
source rows, 19,360,678 <=.5; these are not the study denominator.

Observed shared comparison calls=1, metric aggregation=0, legacy classifiers=0,
AF workers=0, strand workers=0, rules evaluated=False, no approval and no automatic
replacement. Pearson stopped before metrics; no ranks were requested. Exact error:

```text
Whole-study external AF validation stopped: insufficient_overlap_fraction. No chromosome analysis was started; no automatic replacement or dataset retry will be attempted.
```

The old margin reason was absent. All original shard bytes and scientific rows,
all eleven payload hashes and all 403 source inventory/hash entries were unchanged.
Temporary shards/spools were removed. The harness exit was 0 because it verified
the expected stopped analysis; **the external source did not pass harmonisation**.
The direct global-entry harness does not execute the service/retry controller;
separate real service regressions above establish that boundary.

Exact final command (real-data runtime has thread-limit variables unset):

```bash
PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python /Users/JJOHN41/Documents/developing_software/postgwas_v2/.superpowers/sdd/2026-11-05-whole-study-external-af/global-real-data-reproduce.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --source-manifest /Users/JJOHN41/Documents/developing_software/postgwas_v2/.superpowers/sdd/2026-11-05-whole-study-external-af/task-1-fix1-source.json --freeze-id final-fix1-src-8bb78f6f --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f --report /private/tmp/postgwas-whole-external-af-nkifuy6f/global-real-data-final-result.json
```

Harness SHA256 `bec5bd54a251cc85e0011dd3bdfe62fbfd906077b370f899853ee13bcc895ba6`;
final report SHA256 `076fbbc5580f0a21aeb55f33030b5924638f9cb7f5e561abf83d56316b21f091`.
The prior extraction receipt records full-source size/MD5/gzip verification of an
11,749,765-row cloud object generation `1793745884817769`; current cloud/server
identity was not reverified. This local chr6/7/19/20/22 control supports no full22,
all-11-workflow, server-reference-equivalence or successful-acceptance claim.

### Independent generated preparation/adoption comparison

Fresh baseline/final processes used the same immutable two-chromosome corpora at
4,000, 100,000 and 1,000,000 rows. Baseline source SHA256:
`3aa1ec893f4bc18091afca41a05f03fa7a5a194380b7ff7731c96982b59e7016`.
Harness SHA256:
`b30bfd153e62bfdf5b6d64923f7e909543ff008981a615bcdf470181bb55789e`.
Snapshot/live source, harness and corpus fingerprints matched after all runs.
Each baseline completed before its corresponding final, then the comparison ran;
measurements were sequential after Task 1 covering tests completed.

The independent corpus uses no native AF, stable source IDs, permuted positions,
all four study strand orientations, independent external swaps, and 5% sites with
a second ALT AF=.05. Of each 1,000 primary rows, 940 have variable minor AF, 20
have x=.49/r=.51 and 40 have x=.88/r=.90 in external counted-allele direction.
All second-ALT sites are external-swapped; complete-site final REF AF is
`1-(x+.05)`. Production matching and new focused fixtures do not supply the oracle.
Every output row was checked independently for source ID, chr, pos, final alleles,
AF, beta, Z, strand action, reference AF and AF difference; every value/order/dtype
also equals baseline exactly. An initial oracle's `(1-x)-.05` caused 92 exact-double
mismatches; complementing the complete ALT mass `1-(x+.05)` fixed operation-order
parity. The corrected unchanged harness passed positive comparison and changed-AF,
changed usable-row-QC and changed-global-correlation negative self-checks.

Diagnostic pair counts are 3,800/95,000/950,000 (95% fixed-N coverage) with
200/5,000/50,000 known multiallelic exclusions. Independent truth-array direct
r≈.9998210184, MAE=.0012, minor profile=.94 and major profile=.06 agree with final
R1=True/R2=False evidence. All three correlations/errors, raw counts, profiles,
canonical thresholds, exclusions and rule outcomes are separately checked; the
only numeric tolerance is rtol=1e-12/atol=1e-14 for independent floating metrics.
Frames remain exact. Per chromosome at every size, 221 QC leaves are identical,
65 enumerated semantic/accounting changes are allowed, and unexpected changes=0.
No blanket QC exclusion hides scientific output differences.

| Rows | Baseline AF seconds | Final AF seconds | Baseline sampled RSS bytes | Final sampled RSS bytes | Added prepared-partition bytes |
|---|---:|---:|---:|---:|---:|
| 4,000 | 0.143166 | 0.540187 | 145,424,384 | 181,272,576 | 52,211 |
| 100,000 | 0.630234 | 2.295481 | 259,997,696 | 696,614,912 | 471,833 |
| 1,000,000 | 6.058029 | 18.228685 | 987,430,912 | 2,931,163,136 | 4,638,717 |

These are **functional preflight overhead observations, not an optimization**.
AF seconds sum preparation and worker AF phases. Total external parses/joins are
baseline=2/final=2, default parses/orientation matches=2/2. The new path moves them
into preparation, performs one pooled semantic comparison (baseline two legacy
comparisons), and all final workers perform **zero** instrumented reference-table
parses, heavy matches or semantic classifications. This initial-external route
does not have the earlier replacement route's claimed 2-to-1 join benefit.
Added bytes combine default candidates and assignment lists, not isolated AF cache.
RSS was sampled every .02s, may miss peaks, and is not a hard memory limit. Imports,
source hashing, result writes and oracle/comparison time are excluded from measured
AF phases; single runs, fixed baseline-before-final order and uncontrolled OS caches
limit performance interpretation. No 100–200M capacity, full pipeline or server/Linux
claim follows. Any acceptable runtime/RSS decision remains explicit future work.

Exact final nine-command sequence, run from repository root with analysis bin
first on PATH, `PYTHONDONTWRITEBYTECODE=1` and `POLARS_MAX_THREADS=4`:

```bash
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/baseline/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_4000 --output baseline_fix1_generated_4000 --role baseline --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_4000 --output final_fix1_generated_4000_8bb78f6f --role final --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f compare --baseline baseline_fix1_generated_4000 --final final_fix1_generated_4000_8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/baseline/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_100000 --output baseline_fix1_generated_100000 --role baseline --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_100000 --output final_fix1_generated_100000_8bb78f6f --role final --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f compare --baseline baseline_fix1_generated_100000 --final final_fix1_generated_100000_8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/baseline/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_1000000 --output baseline_fix1_generated_1000000 --role baseline --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f run --corpus generated_1000000 --output final_fix1_generated_1000000_8bb78f6f --role final --freeze-id formal-round1-src-8bb78f6f
/Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B .superpowers/sdd/2026-11-05-whole-study-external-af/generated-af-compare.py --source-root /private/tmp/postgwas-whole-external-af-nkifuy6f/final-fix1-8bb78f6f/src --work-root /private/tmp/postgwas-whole-external-af-nkifuy6f compare --baseline baseline_fix1_generated_1000000 --final final_fix1_generated_1000000_8bb78f6f
```

### Root regression/static checks and explicit baseline limitations

| Check | Observed final result |
|---|---|
| Broad workspace, exact 34 owners below | 1,185 passed, 159 warnings, 460.88s; exit 0 |
| Clean Git-base archive plus only 14 owned paths, exact 10 portable owners | 492 passed, 266.01s, no warnings; exit 0 |
| Unrestricted pytest | 37 existing collection errors, one warning, 10.94s; exit 2 |
| Exact unrestricted error-set comparison | All 37 match earlier permanent baseline; no new/resolved paths |
| AST/configuration/protected-path/whitespace checks | Eight changed/new Python owners parse; all non-help/origin YAML identical to base; protected adapters unchanged; diff check passes |
| Generated policy / wiki | Policy synchronization passes; 50-page wiki check passes |
| Wiki CLI example extractor | Exit 1, exactly two known installer false positives below |

The 159 broad warnings originate at nine byte-identical-to-base producer locations:
seven concordance joins in `concordance/analysis.py:205,814,822,1071,1096,1624,1641`,
`summary_statistics_io.py:208` and `test_harmonisation_metadata_preamble.py:178`.
No warning producer was edited. This source comparison does not claim identical
warning counts to a differently scoped baseline or that the warnings are fixed.
Unrestricted's one warning is the third-party marshmallow distutils deprecation.

The complete exact 37-path collection-error list remains in the permanent
[earlier baseline record](2026-11-04-replacement-af-validation.md#explicit-unresolved-baseline-checks).
Current `unrestricted-baseline-comparison.json` confirms exact set equality, not
merely equal counts. Failures include archived resource-test copies, the current
resource-preparation tests and `tests/test_harmonisation_coordinates.py`.
The CLI checker currently reports only `README.md` and
`docs/wiki/getting-started/installation.md`: unsupported command `--jobs` from
`postgwas --jobs 2`. Its extractor mistakes valid installer `--name postgwas --jobs 2`
arguments for application commands. Those sources and unrelated tooling were not
changed. These issues remain disclosed, so no repository-wide release-green claim
is made. Local broad owners can use existing ignored helpers/private history;
the isolated candidate is the clean affected-portable-owner check, not a full fresh
clone or all-test proof. No unrelated collection/CLI repair, ignore-policy change
or raw scientific fixture is included.

Exact root commands below use
`PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH`,
`PYTHONDONTWRITEBYTECODE=1`, `POLARS_MAX_THREADS=4`, and the cwd-specific absolute
`PYTHONPATH` shown. The broad/unrestricted cwd is this repository; portable cwd is
`/private/tmp/postgwas-whole-external-af-nkifuy6f/verification/candidate` (Git archive
of base plus exactly the 14 owned paths). The three root suites ran concurrently
after controlled generated measurements finished; their elapsed times are not
benchmarks. Full outputs and argv were inspected, not rerun for this prose record. Root also
verified that all nine owned runtime/configuration/test files are byte-identical
to the tested clean candidate and all 403 source files still match the final freeze.

```bash
# broad; cwd: /Users/JJOHN41/Documents/developing_software/postgwas_v2
PYTHONPATH=/Users/JJOHN41/Documents/developing_software/postgwas_v2/src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_delimiter_header_selection.py tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_build_fallback.py tests/test_harmonisation_config.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_exclude_chr_x.py tests/test_harmonisation_exclude_x_accounting.py tests/test_harmonisation_exclude_x_interface.py tests/test_harmonisation_exclude_x_references.py tests/test_harmonisation_exponentiated_se.py tests/test_harmonisation_exponentiated_se_integration.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_metadata_preamble.py tests/test_harmonisation_multiallelic_reference_frequency.py tests/test_harmonisation_reference_frequency_export.py tests/test_harmonisation_reference_frequency_provenance.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_resource_preflight.py tests/test_harmonisation_sample_sheet_generator.py tests/test_harmonisation_shared_coordinate_identifier.py tests/test_harmonisation_study_report_integration.py tests/test_harmonisation_allele_join.py tests/test_harmonisation_imputation_quality.py tests/test_harmonisation_eaf_replacement_downstream.py tests/test_harmonisation_chromosome_partition_index.py tests/test_harmonisation_strand.py tests/test_harmonisation_external_reference_staging.py tests/test_harmonisation_adaptive_scheduling.py tests/test_harmonisation_external_af_preflight.py tests/test_harmonisation_allele_frequency.py --basetemp /private/tmp/postgwas-whole-external-af-nkifuy6f/verification/fixtures-broad
```

```bash
# portable; cwd: /private/tmp/postgwas-whole-external-af-nkifuy6f/verification/candidate
PYTHONPATH=/private/tmp/postgwas-whole-external-af-nkifuy6f/verification/candidate/src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_resource_preflight.py tests/test_harmonisation_external_af_preflight.py --basetemp /private/tmp/postgwas-whole-external-af-nkifuy6f/verification/fixtures-portable
```

```bash
# unrestricted; cwd: /Users/JJOHN41/Documents/developing_software/postgwas_v2
PYTHONPATH=/Users/JJOHN41/Documents/developing_software/postgwas_v2/src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider --basetemp /private/tmp/postgwas-whole-external-af-nkifuy6f/verification/fixtures-unrestricted
```

Static/document command owners:

```bash
PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B tools/docs/update_harmonisation_policies.py --check
PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B tools/docs/build_wiki.py --check --output /private/tmp/postgwas-whole-external-af-nkifuy6f/wiki-check
PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B tools/docs/validate_wiki_cli.py
git diff --check
```

The AST owners are the seven changed Python source files plus the focused new
test listed in scope. Recursive YAML comparison omits only mapping keys named
`help` and `origin`. Source/focused-test hashes and all non-document owned hashes
were rechecked after this prose change; whitespace and relative-document links
were checked without repeating scientific suites. No separate repository-wide
linter/typechecker run is claimed.

### All seven ledger rulings and their costs

1. Ruling: Reuse the current checkout on a new codex feature branch, preserving unrelated files — existing workspace is suitable — wrong attribution costs restaging/rework.
2. Ruling: The user's correction request authorizes extending the existing whole-file design to the missing-study-AF path without another approval pause — this is the previously omitted source mode — wrong scope costs visible reviewable rework.
3. Ruling: Pool the raw source screen over observed study chromosomes after existing deduplication/staging, with counts rather than averaged fractions — this matches dataset scope and bounded existing I/O — wrong scope could affect screening for partial chromosome datasets; explicit scope tests/documentation required.
4. Ruling: Preserve the approved fixed study-row denominator and all thresholds even though the real five-chromosome fixture has only .46818386144455965 eligible coverage — do not manufacture acceptance — users may still need a suitable source or explicit future coverage-policy decision; full22 success is not asserted.
5. Ruling: Rule2-only external MAF support stops without selecting default or explicit replacement — original study AF needed for the native replacement contract is absent — a different desired recovery policy would require a separately specified change.
6. Ruling: Keep direct helper legacy compatibility but require the new bound approval in normal dataset workers — no public helper regression while the real workflow cannot fall back — missing integration coverage could leak the old classifier, requiring a real service boundary regression.
7. Ruling: Existing unrestricted collection/CLI-extractor issues stay disclosed and outside scope; fix every in-scope regression — do not alter user artifacts or unrelated tooling — no repository-wide release-green claim.

### Compliance and remaining work at this snapshot

Checklist: **passed for scoped implementation, correction review and the evidence
record; final whole-task review/commit/cleanup pending**. Cumulative diff reviewed:
**yes**, all 14 owned paths audited against base/current freezes and the complete
original/correction review reports; the two-document prose delta checked separately.
Hardcoding and duplication audit: **passed** for scoped code/docs; no new parameter,
resource or threshold/default owner. Scientific validation: **completed for the
recorded synthetic, compatibility and local insufficient-coverage boundaries**.
Tests run: exact commands, data, expected/observed results and owners above.
Tests not run: full22/all-11 workflows, native complete PostGWAS production pipeline,
100–200M rows, production/Linux server, renewed cloud/server provenance, repeated
timing trials and repository-wide lint/type checks; these lack execution evidence
and remain unverified. Outstanding issues: known 37 collection errors/two extractor
false positives, measured preflight overhead and the intentionally below-gate real
fixture; none is silently promoted to success. Final whole-task review, authorized
exact-path local commit and owned scratch cleanup remain root-owned. User `_to_delete/`,
`tmp/`, duplicate ` 2` copies and unrelated ignored fixtures remain preserved.
No deferred minor or additional ruling is introduced by this documentation cycle.

Evidence provenance before cleanup: `task-1-report.md`, `task-1-fix-report.md`,
`task-1-review.md`, `task-1-fix1-review.md`, `task-1-fix1-source.json`,
`generated-final-{summary,commands,source}.json`, `generated-final-{4000,100000,1000000}.json`,
`generated-validation.md`, `real-data-baseline.json`, `real-data-recon.md`,
`global-real-data-{commands,final-result}.json`, `global-real-data-final-validation.md`,
`root-verification-{commands,results}.json`, `post-test-source-verification.json`,
broad/portable/unrestricted logs,
`static-results.json`, `warning-source-comparison.json`,
`unrestricted-baseline-comparison.json`, and `wiki-cli-check.log` in the private task
workspace. Material outcomes and commands are preserved here before its removal.

Accuracy confidence: **98%** for this verified pre-commit record; broad and portable
results, current source identity and independent row/numeric oracles agree, while
known collection/tooling defects, single-run resource measurements and unperformed
whole-workflow/server/capacity checks bound the claim.

### Final wording correction gate — 2026-11-05

The independent whole-task review found only two P3 wording issues. In the shared
pair-count guard, `Replacement AF replay` became source-neutral `AF replay`.
The private `_prepare_only` docstring now names both initial-external and explicit-
replacement preparation. The exception type and trigger, numeric calculations,
configuration, defaults, source selection and publication controls are unchanged.

The measured scientific and regression evidence above belongs to the exact
`8bb78f6f31e65f8264baf1c0811c590a9b7d25f104b7be3e71e58ba288aafdfa`
source freeze. The 1,185-pass broad, 492-pass portable, 284-pass covering,
five-chromosome and generated 4K/100K/1M results were not rerun on the final
wording bytes. A focused final gate exercised seven collected cases and passed;
the exact command and AST proof are in the private `final-fix-report.md`.

Final wording source: **403 regular `src` files**, SHA256
`e06b266b7e06cccf14aebba576660ab12bba4a37d34935677d6feb6d5b54e6c5`;
source plus the focused test SHA256
`c7024ab5e25c5a8f00baa919850bd6797e29b3b25b6db3d32edc2bc16d78e1b5`.
The private `final-fix-source.json` records every per-file identity. Its ordering
is the same component-wise `Path` order as the measured freeze. Comparing all 403
per-file identities with the measured manifest changes only
`src/postgwas/modules/harmonisation/allele_frequency.py`; normalizing exactly
the two changed strings makes the full Python AST identical. Final scoped review,
exact-path local commit and owned scratch cleanup remain root-owned.
