# PostGWAS `sumstat_filter` — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/filtering/` (`sumstat_filter.py` 1243, `service.py` 179, `cli.py` 104), `config/models/modules/filtering.py`, `config/defaults/modules/filtering.yaml`, `docs/wiki/modules/filtering.md`, `tests/test_filtering_summary.py`, and the module's contract with `qc_summary` and the pipeline.
**Method:** read-only on the tree — no file was modified. Every claim about bcftools behaviour below was **executed** against bcftools 1.19 / htslib 1.19 on purpose-built fixtures, including a side-by-side run of the module's real filter pipeline against its reason-attribution pipeline using the shipped defaults. The polars numbers in **I3** are measured on a 5,000,000-row tag file. Line numbers refer to the working tree as of 2026-08-26.

**Summary.** This is the most carefully engineered module in the review series. It already does the things the previous two reviews had to ask for: it validates the VCF header contract with the required fields derived from the *active* rules, infers the build only from `##genome_build`, matches MHC contig naming against the header and errors when it cannot, converts BED coordinates correctly, runs `set -euo pipefail` with a post-hoc index check, keeps `--threads` off the uncompressed intermediates, and refuses to print a fake `0` for a count it could not obtain. Most of what follows is about the *audit* the module runs alongside the filter — which uses a different operator than the filter itself, costs a second full pass over the VCF, and, when it disagrees with reality, says so and ships anyway.

---

## CRITICAL

### C1 — The reason audit and the real filter use two different non-SNP operators; the module detects the resulting mismatch and publishes anyway

**Where:** `sumstat_filter.py:977-978` (`--types snps`) vs `:896-901` (`TYPE != 'snp'`); reconciliation at `:1095-1109`.

The real pipeline removes non-SNPs with a bcftools *option*:

```python
cmd1 = "%s view -Ou" % quoted_bcftools
if include_expr is not None: cmd1 += " -i %s" % shlex.quote(include_expr)
if not include_indels:       cmd1 += " --types snps"
```

The audit attributes that same removal with a filter *expression*:

```python
condition_checks.append({"label": "Indels and other non-SNP variants",
                         "expr": f"({variant_type_field} != 'snp')", "priority": 50})
```

**These are not the same operator.** Verified on bcftools 1.19 with a record `A > T,AT`:

| record | `view --types snps` | `-i "TYPE='snp'"` | tagged by `-e "TYPE!='snp'"` |
|---|---|---|---|
| `A > G` (SNP) | kept | kept | no |
| `A > AT` (insertion) | dropped | dropped | yes |
| `A > T,AT` (mixed multi-allelic) | **kept** | **dropped** | **yes** |

`--types snps` keeps a site if *any* alternate allele is a SNP; the expression tags it as a non-SNP removal. Running the module's full default pipeline against its full default audit on a 12-record fixture:

```
real pipeline   : before=12  after=4   actual_removed=8   (the multi-allelic survives)
reason pipeline : primary_removed_total=9                 (the multi-allelic is attributed to "Indels and other non-SNP variants")
```

**The module notices and continues.** `reason_statistics["reconciled"]` goes False (`:1097-1100`), a warning is logged, `Reason reconciliation … FAILED` is written into the TSV (`:361-362`), the screen prints a red `Removal count check: 9 assigned to reasons != 8 removed in total` (`:166-175`) — and the filtered VCF and the report are then published unchanged (`:1136-1146`). Compare `qc_summary`, whose analogous accounting check *raises* `VcfAssessmentError` and refuses the report (`qc_summary/assessment.py:570-581`). Two sibling modules, the same invariant, opposite consequences.

**Scope, honestly.** Harmonisation splits multi-allelics with `bcftools norm -m-any` (`harmonisation/vcf_processing.py:495, 629`), so for a VCF produced by this pipeline the mismatch is latent. It is live for the externally supplied VCFs the standalone `postgwas sumstat_filter` command exists to serve. A second, guaranteed reconciliation failure needs no multi-allelics at all: `empty_expression_action: match_none` with no active rule sets `include_expr = "1"`, which — verified — matches **nothing** (12 in, 0 out, exit 0), removing every variant with zero attributed reasons.

**Fix.** Use one operator in both places. Replacing `--types snps` in stage 1 with `-e "TYPE!='snp'"` makes the attribution exact by construction and costs nothing. Then make an unreconciled audit fail the step: a report the module already knows is wrong should not be written as a result.

### C2 — On a multi-sample VCF every FORMAT-field rule is satisfied by whichever study is most permissive, and nothing checks the sample count

**Where:** `sumstat_filter.py:711-807` (all `FORMAT/…` expressions); no `bcftools query -l` and no `--samples` anywhere in `modules/filtering/`.

Every threshold that matters — MAF, imputation quality, significance — is expressed against a FORMAT field. On a site-level `-i`, bcftools keeps the record if **any** sample satisfies it. Verified on a two-sample VCF:

| record | STUDY_A AF | STUDY_B AF | kept by `-i 'FMT/AF>=0.01 & FMT/AF<=0.99'` |
|---|---|---|---|
| `both_ok` | 0.30 | 0.30 | kept |
| `A_fails_maf` | **0.001** | 0.30 | **kept** |
| `B_fails_maf` | 0.30 | **0.001** | **kept** |
| `both_fail` | 0.001 | 0.002 | dropped |

So a variant that fails MAF in the study of interest survives because some other study in the file passes. The output VCF is then trusted by every downstream module.

**The obvious fix does not work.** Adding `-s` to the same `view` still keeps `A_fails_maf` — bcftools evaluates the expression before subsetting within one invocation. Verified that the subset must be its own stage:

```
bcftools view -H -s STUDY_A -i 'FMT/AF>=…' multi.vcf.gz     → both_ok A_fails_maf B_fails_maf A_fails_info   (wrong)
bcftools view -Ou -s STUDY_A … | bcftools view -H -i 'FMT/AF>=…'  → both_ok B_fails_maf A_fails_info          (correct)
```

`qc_summary` handles this properly: it resolves the study sample with `select_vcf_sample` and passes `--samples` to its query. The module that only *counts* is careful about sample identity; the module that *writes the VCF* is not.

The docs do list this under Limitations ("does not select one study/sample from a multi-study GWAS-VCF"), so it is known — but the wording understates it, and disclosure does not stop a silently wrong VCF. A one-command guard (`bcftools query -l | wc -l`, refuse or require an explicit sample) converts a silent wrong answer into an error; full support is the separate-stage subset above.

---

## IMPORTANT

### I1 — `info_max` cannot be set to the value the rest of the codebase expects

**Where:** `config/models/modules/filtering.py:156` — `info_max: float | None = Field(default=None, ge=0, le=1)`; compare `config/models/modules/qc_summary.py:45` — `info_max: float = Field(ge=0, le=100)` with the shipped default `1.05`.

`modules.filtering.info_max: 1.05` **fails schema validation**. Two consequences:

- The QC preview reports `Imputation quality outside 0.7 ≤ SI ≤ 1.05` and the filter cannot be configured to match it, so the preview and the filter necessarily disagree on the upper bound.
- With the shipped `info_max: null` there is no upper bound at all, so an SI of 3.0 — which does not mean "very well imputed", it means the field is corrupt or mis-mapped — passes the filter while failing the QC preview.

Minimac/MACH r² legitimately reaches slightly above 1, which is exactly why harmonisation carries a purpose-built warning comparing `info.mach_rsq_max` against `modules.qc_summary.rules.info_max` (`harmonisation/service.py:4660-4670`). Raise filtering's ceiling to match its sibling and give it a real default.

### I2 — The audit is a second full pass through twelve chained bcftools processes, and its failure path costs ten more

**Where:** `sumstat_filter.py:276-298`; fallback at `:943-948`.

```python
stages = ["bcftools view -Ou <input>"]
for check in checks:
    stages.append("bcftools filter -Ou -m + -s <TAG> -e <EXPR>")
stages.append("bcftools query -f '%FILTER\\n'")
```

Under the shipped defaults there are ten checks (missing `FORMAT/AF`; MAF; missing `FORMAT/SI`; INFO range; missing `INFO/AF`; missing `INFO/EUR`; AF difference; non-SNP; palindromic; MHC), so this is one `view`, ten `filter` and one `query` — twelve processes, each decoding and re-encoding BCF over the entire file — run **in addition to** the real filtering pipeline. The module therefore reads the VCF end to end twice.

When that pass fails, the fallback calls `count_matching` once **per check** — each a full `bcftools query -i EXPR | wc -l` scan — so the degraded path is ten more complete passes, and it is reached through a `log_warn` that a user reading the terminal will not notice.

The tagged stream already contains everything the output needs: "removed" is exactly "carries at least one removing tag". One pass can `tee` the tagged BCF into (a) `query -f '%FILTER\n'` for the audit and (b) `view -e '<union of removing tags>' -Oz` for the filtered VCF. That halves the I/O **and** makes reconciliation true by construction — which is also C1's fix.

### I3 — The polars aggregation is neither streaming nor cheap: measured 2.8× slower and 3.9× more memory than necessary

**Where:** `sumstat_filter.py:191-229`.

```python
scan = pl.scan_csv(tag_file, has_header=False, separator="\t", schema={"filter_tags": pl.String})
tags  = pl.col("filter_tags").fill_null("").str.split(";")
masks = dict((check["tag"], tags.list.contains(check["tag"])) for check in checks)
...
values = scan.select(selections).collect()        # no engine="streaming"
```

`str.split(";")` materialises a `List[String]` column of table height before ten `list.contains` probes run over it, and `.collect()` without `engine="streaming"` holds it. Measured on 5,000,000 rows / 57 MB with ten tags:

| implementation | wall time | peak RSS |
|---|---|---|
| current — `str.split` + `list.contains`, `.collect()` | **16.3 s** | **555 MB** |
| `str.split` + `list.contains`, `collect(engine="streaming")` | 16.9 s | 146 MB |
| `str.contains(literal=True)`, `.collect()` | 5.2 s | 344 MB |
| `str.contains(literal=True)`, `collect(engine="streaming")` | **5.9 s** | **142 MB** |

All four returned identical counts (`primary_removed_total`, `overlap_variants`, `extra_rule_matches` and every per-tag figure). `engine="streaming"` alone is a 3.8× memory cut; the literal substring match is a 2.8× time cut; together it is a two-line change.

One caveat to record in a comment if this is taken: `str.contains` is a substring test, so it is exact only while no tag name is a prefix of another. The `%02d` tag format (`sumstat_filter.py:265`) guarantees that up to 99 checks; a hypothetical hundredth would make `PGWAS_FILTER_10` a prefix of `PGWAS_FILTER_100`.

### I4 — `--threads` never reaches polars

`execution.threads` is threaded into the bcftools stages correctly — and, to the module's credit, correctly *omitted* from the `-Ou` intermediates and applied only to the final `-Oz` (`:973, 983, 1011, 1030-1031`). But `_summarize_filter_tags` sizes its thread pool from the machine's core count. `POLARS_MAX_THREADS` appears exactly once in `src/`, in `harmonisation/service.py:224`, with a comment explaining it must be set before the process starts. Same finding as `qc_summary` I4, same one-line fix.

### I5 — No completion manifest, no resume, no resolved-configuration file, and the log is written only at the end

`run_sumstat_filter_direct` has none of the `resolve_completion_resume` / `write_completion_manifest` / `write_resolved_configuration` machinery that `qc_summary` uses, and the `sumstat_filter` registry entry carries no `direct_checkpoint` (`qc_summary`'s is `"native"`, `registry.py:428`). So the most expensive module in this group is the one that cannot be resumed.

It also builds its own logging from a `StringIO` plus `log_print`/`log_warn` closures (`:497-526`) rather than the shared `PipelineLogger`, and flushes it exactly once, at the end, via `write_log()`. The error paths call `write_log()` too — but a hard kill during the long bcftools pipeline (OOM, cluster preemption, node failure) loses the entire log, including the resolved configuration and the exact command line, which are the two things needed to diagnose it.

---

## OPTIONAL

**O1 — an entire report section is unreachable.** `_filtering_summary_lines` renders an eight-field "Variant flow" block from `data_flow` (`:74-91`), and the only production call site (`:1127`) never passes `data_flow`. Its sole exercise is `tests/test_filtering_summary.py:65`. Either wire the harmonisation pre-VCF counts through or delete the parameter — right now a test is keeping dead code alive.

**O2 — "unavailable" is reported as "FAILED".** When the tagging pass fails, `report_statistics` is constructed with `"reconciled": False` (`:1116`), so the TSV records `Reason reconciliation … FAILED` even though nothing was mis-attributed — it simply was never measured. This module is otherwise scrupulous about that distinction; `_fmt_count`'s docstring is literally "never a fake 0". Use a third state.

**O3 — two implementations of the same report primitives.** `_is_missing`/`_fmt_count` (`:30-47`) duplicate `metric_available`/`format_metric_count`, which `qc_summary/reporting.py` already exports in `__all__`. `_filtering_summary_lines` likewise re-implements the flow / conditions / accounting / reconciliation layout that `qc_summary_lines` builds — including the same harmonisation `total_variant_*` keys.

**O4 — a temporary file is created only to mint a path, deleted, then written twice.** `:1118-1126` opens a `NamedTemporaryFile`, closes it and immediately `unlink()`s it; `_write_filter_reason_report` then creates its *own* temporary file and `os.replace`s onto that path; `:1145` replaces it again onto the destination. One atomic write is enough.

**O5 — two of the three missing-value actions have no CLI flag.** `--missing-info-action` exists (`common.py:1233`) and is mapped in `MODULE_CLI_OVERRIDES`; `missing_af_action` and `missing_pvalue_action` are YAML-only and absent from the map. They are equally scientific policies.

**O6 — `missing_pvalue_action` is inert under the shipped defaults.** The LP missing check is appended only when `minimum_neglog10_p is not None` (`:711-719`), and the default is `null`, so `missing_pvalue_action: remove` can never fire. The behaviour is right; the default reads as an active policy.

**O7 — pre- and post-imputation filtering share one configuration.** `post_imputation_filter` reuses `run_sumstat_filter_runner` against the same `modules.filtering` block (`registry.py:146-158`), so both passes apply identical MAF, palindromic, MHC and variant-type rules. Only the INFO rule is scientifically meaningful the second time — imputation changes SI — while the rest re-remove nothing and still pay the full double-pass cost of I2. There is no way to configure the two stages differently.

**O8 — the summary always prints to stdout.** `print(line)` at `:1174-1176`, unconditionally. `qc_summary` gates the equivalent on `configuration.logging.show_screen`.

**O9 — one call uses the unresolved input path.** `count_matching(vcf_path, check["expr"])` at `:947` passes the raw parameter, while every other use is `str(input_vcf)`, the expanded and resolved path.

---

## Verified correct — do not re-flag

All executed against bcftools 1.19 / htslib 1.19 with purpose-built fixtures.

- **`-i '1'` matches nothing** — 12 records in, 0 out, exit 0, empty stderr. The comment at `:844-846` is right, and omitting `-i` entirely really is the only way to express "match everything". The `empty_expression_action` handling and its warning are correct.
- **`FMT/X == '.'` correctly identifies a missing FORMAT value** with the single-quoted literal exactly as written, and `!= '.'` is its exact complement. Every `keep`/`remove` combination in `:711-842` was traced on the fixture and the resulting tags matched the include expressions.
- **`-T ^bed` silently excludes nothing on a contig-name mismatch** — a `chr6` BED against a `6` VCF kept 12/12 with exit 0 and empty stderr. `_match_contig_naming` (`:1206-1243`) is a real guard against a real, silent failure; it accepts the equivalent spelling with a warning and *errors* when neither matches. This is precisely the guard `annot_ldblock` lacks.
- **The BED start decrement is necessary and correct.** Excluding 1-based position 30,000,000 requires BED start 29,999,999; with the un-decremented start the first base of the region survives. Verified in both directions. The comment at `:989-991` explains exactly the bug it prevents.
- **Thread placement is right** — no `--threads` on the `-Ou` intermediate stages, threads only on the final `-Oz` compression.
- **`set -euo pipefail` on both shell pipelines**, plus a post-hoc check that the temporary VCF *and* its index exist and are non-empty (`:1044-1065`). This is what stops a mid-pipe abort from being reported as "variants removed by filtering", and the comment says so.
- **Build inference from `##genome_build` only** — exactly one declaration required, filename / `##reference` / contig-`assembly` explicitly rejected — with `required_fields` derived from the *active* rules so an unused tag is not demanded (`:544-558`). This is the best input validation in the modules reviewed so far and is what `qc_summary` and `annot_ldblock` should copy.
- **MHC intervals** match the GRC regions for both builds and match `qc_summary`'s.
- **Atomic finalisation** — VCF, index, MHC BED and reason report are all staged and `os.replace`d, with cleanup on every failure branch (`:1136-1159`).
- **First-failure attribution is deterministic and visible**: checks carry explicit `priority` values, are sorted before both display and attribution (`:916-932`), and the report states the order. The `seen`/`primary` accumulation in `_summarize_filter_tags` implements it correctly — confirmed against the module's own unit test and re-derived on the 5M-row benchmark.

---

## Documentation

`docs/wiki/modules/filtering.md` is accurate and unusually precise — the `##genome_build` policy, the MHC BED conversion, and the `7.30103` LP-versus-P note are all correct and worth keeping. Three gaps:

- The Limitations note about multi-study VCFs understates C2. The behaviour is not "does not select one study/sample"; it is "every FORMAT-based rule is satisfied by whichever study in the file is most permissive".
- Nothing documents what happens when reconciliation fails — that the filtered VCF is still published and the per-reason counts are known to be wrong.
- Common problems says missing values "are never silently treated as passing values". That holds for the three `missing_*_action` rules but not for the palindromic exclusion: `-e (palindromic & AF >= lower & AF <= upper)` is false when `FORMAT/AF` is missing, so a strand-ambiguous SNP with no frequency is always kept, whatever `missing_af_action` says. Reachable whenever `maf_min` is null or `missing_af_action: keep` — and the same gap exists in `qc_summary`'s palindromic rule.

---

## Suggested order of work

1. **C1** — use one non-SNP operator in both places (`-e "TYPE!='snp'"` in stage 1), and make an unreconciled audit fail the step rather than publish. Small change; removes a class of wrong report.
2. **C2** — a `bcftools query -l` guard that refuses a multi-sample VCF, or a separate `view -s` stage (verified: it must be its own stage).
3. **I1** — raise `filtering.info_max` to match `qc_summary` and give it a default, so the preview and the filter can agree.
4. **I3** — `engine="streaming"` plus `str.contains(literal=True)`. Measured 2.8× faster, 3.9× less peak memory, identical numbers, two lines.
5. **I4** — `POLARS_MAX_THREADS`, copying `harmonisation/service.py:221-224`.
6. **I2** — collapse the audit and the filter into one tagged stream. The largest change here; it subsumes C1's fix and halves the module's I/O.
7. **I5** — completion manifest, resume, and the shared `PipelineLogger`.
8. **O1–O9.**

Steps 1–3 fix numbers that are wrong. Steps 4–6 are cost. Step 7 is the operational gap that will be felt on the first preempted job.
