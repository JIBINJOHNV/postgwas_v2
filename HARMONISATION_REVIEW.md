# PostGWAS harmonisation module — code and documentation review

Reviewed: `src/postgwas/modules/harmonisation/` (53 files, ~33k lines), its config surface
(`config/models/modules/harmonisation.py`, `config/defaults/modules/harmonisation.yaml`),
`README.md`, and `docs/wiki/harmonisation/*`.
Every finding below was read in the source and, for the Critical items, adversarially
re-verified against alternative code paths and config defaults.

**Overall verdict.** The step ordering is scientifically correct and the module is unusually
well instrumented — per-variant reject provenance, a per-chromosome row-count reconciliation
assert, `DECIDE` logging of study-level inferences, and tail-safe p↔Z conversions via
`log_ndtr`/`ndtri_exp`. The classic sign traps are handled correctly (see
[Verified correct](#verified-correct)). The problems are concentrated in three places:

1. **Shipped policy defaults that are permissive where the check matters most** — four of them
   let scientifically wrong values reach the released VCF with no hard stop.
2. **Two study-level inferences that have no data-based cross-check** (EAF orientation,
   binary-trait Neff), even though the evidence needed to check them is already loaded.
3. **The dataset stage is fully eager Polars** while the `concordance/` sub-module already
   demonstrates the correct streaming pattern. This is the scalability ceiling.

Nothing below asks for a refactor for its own sake. Where consolidation is proposed, it is
because two copies have already diverged in behaviour.

---

## Actual workflow order (ground truth)

This is what the code runs, and it is the order the documentation should describe.

**Dataset stage** — `run_harmonisation_pipeline` (`service.py:3810`) → `_prepare_dataset_for_chromosome_processing` (`service.py:2144`)

| # | Step | Scientific purpose |
|---|---|---|
| 1 | Check the run configuration | Refuse an unusable column mapping before reading a large file |
| 2 | Check the file header | Every configured column exists |
| — | Truncation + data-line count | Unnumbered; runs after step 2's log context closes (`service.py:4110`) |
| 3 | Read the summary statistics | Parse → read-mandatory gate → coordinate/allele normalisation and rejection → numeric cast → **second** coordinate/allele rejection pass → duplicate resolution → internal INFO normalisation |
| 4 | Check the data behind the configuration | `prepare_missing_sample_sizes`, `validate_content` |
| 5 | Genome build | Decided from the data, not from the declaration (`genome_build.py:187`) |
| 6 | Strand detection | Whole-study forward/reverse consensus from non-palindromic reference matches |
| 7 | Study properties | Effect type, SE scale, p-value scale, EAF-vs-MAF |
| — | Resource preflight | Unnumbered; exact per-chromosome resource check (`service.py:2540`) |
| 8 | Split by chromosome + stage external EAF/INFO | |

**Per-chromosome stage** — `process_one_chromosome` (`service.py:1436`)

`01` load partition · `02` reference resources · `03` effect type → log-odds BETA · `04` allele
orientation + EAF · `05` sample size / Neff · `06` BETA+SE from Z · `07` p-value scale ·
`08` SE from BETA+P · `09` Z = BETA/SE · `10` effect-statistics validation · `11` INFO ·
`12` SNP identifier · `13` final completeness gate · *(unnumbered row reconciliation,
`service.py:1849`)* · `14` export adapter TSV · `15` gwas2vcf · `16` normalise/annotate/liftover

**Post-merge stage** — `_run_post_merge_stages` (`service.py:3391`)

`1` concatenate chromosome VCFs · `2` population-frequency comparison · `3` side-file merge +
cleanup · `4` raw-VCF QC assessment · `5` final QC report and logs

**Ordering assessment: correct.** Effect-type normalisation (03) precedes allele orientation
(04) so the orientation flip acts on a log-odds BETA; orientation precedes Neff (05) and
Z-based reconstruction (06) so `2·EAF·(1−EAF)` uses the aligned frequency; p-value handling
(07) precedes SE-from-p (08) which precedes Z (09) which precedes validation (10). INFO (11)
and identifiers (12) correctly come last, after alleles are final. The one genuine ordering
defect is inside dataset step 3 — see [I6](#i6).

---

## Critical

### C1 — A frequency column reporting the *non-effect* allele is never detected; every exported EAF is inverted

`study_properties.py:753-755` decides `eaf_is_maf = (fraction ≤ 0.5) > 0.95`
(`eaf.maf_decision_cutoff: 0.95`). A column holding the frequency of the **other** allele — a
very common real-world mapping error — is U-shaped, so roughly 50 % of values fall below 0.5
and the decision is `False`. `allele_frequency.py:1222` calls `confirm_eaf_with_reference`
**only** under `if study_says_maf:`, and that line is the function's sole call site. The
`False` path falls through to `qc_info["decision_source"] = "internal_eaf_study_decision"`
(`:1288`) with no reference comparison at all. `strand.py:801-803` then flips the
already-wrong value on swapped rows, so the exported AF is `1−AF` for essentially every
variant.

Nothing stops this. The per-variant `strand.af_discordance_action` gate
(`allele_frequency.py:1415-1557`) defaults to **`warn`** (yaml:2444) at tolerance 0.20, so
common variants near 0.5 are not even counted. `--validate` concordance compares the input to the
output VCF, and for `eaf_is_maf is False` it derives its expectation *from the input's own AF*
(`concordance/analysis.py:1599-1602`: `when(orientation_factor > 0).then(input_af).otherwise(1 - input_af)`).
An inverted column is therefore perfectly self-consistent and PASSes.
`population_frequency.py:45-192` *does* run a dataset-level correlation against AFR/EAS/EUR/SAS
and is enabled by default, but an inverted AF gives r ≈ −1, which falls below
`minimum_correlation: 0.80` (`:147`), so `closest_population` is `None`, the run is reported as
`status: "inconclusive"` (`:179`), and `check_selected_population` maps that to
`"not_compared"` with **no warning at all** — it warns only on `"mismatch"`
(`population_frequency.py:218-230`).

**Why it matters.** GWAS-SSF and GWAS-VCF both require AF aligned to the effect allele. β/SE
reconstruction is unaffected because `2p(1−p)` is symmetric — which is exactly why this stays
invisible — but every downstream MR, PRS and colocalisation allele check is corrupted.

**Change.** Three lines of defence, all cheap because the reference join already happens:
1. Call `confirm_eaf_with_reference` for **every** study frequency column, not only MAF-like
   ones. Report aligned-vs-folded error either way.
2. Add a dataset-level gate on the *fraction* of non-palindromic variants that are AF-discordant
   (fail, or auto-invert with explicit provenance) rather than only a per-variant warn.
3. In `population_frequency.py`, treat a strong **negative** correlation as the distinct,
   named diagnosis "the frequency column appears to be inverted", not as `inconclusive`.

### C2 — Neff is silently set to total N for a case-control study supplying one count column

`sample_size.py:578-593`: with `sample_size.controls_only: use_ncontrol_as_total` (default,
yaml:3668), a study that configures only `ncontrol_col` gets `Neff := that count`. The binary
guard at `:556-563` fires only when `sample_size.trait_type == "binary"`, and the default is
`auto` — resolved **solely** from the sample sheet's declared `trait_type`
(`service.py:835-887`); nothing infers it from the data. The sample-sheet auto-generator
hard-codes `trait_type: auto` (`config/models/modules/harmonisation.py:179`,
`Literal["auto"]`), so generated sheets always land here. There is no `n_col`/`neff_col`
mapping at all (`input_validation.py:76-91`), so a study with a per-variant total-N column
*must* map it to `ncontrol_col`.

`sample_size.py` contains zero references to `effect_type`, even though step 03 has already
resolved `odds_ratio` and written it to `sample_column_dict["effect_type"]`
(`effect_type.py:847`). Nothing downstream catches the result: `effect_validation.py` never
mentions Neff, and `effective_sample_size_above_outlier_threshold` is *read*
(`cli.py:721`, `qc_reporting.py:607`) but never produced anywhere in the tree.

**Why it matters.** Neff is overstated by 2× (balanced design) to 4× (extreme ratio). It feeds
`2·EAF·(1−EAF)·Neff` at `effect_from_z.py:699`, so reconstructed SEs are too small by
√2–2×, and it is written into the output VCF as the N field that LDSC and MiXeR consume for
h² and rg. This is the precise failure the `cases_only` error text at `sample_size.py:596-601`
warns about, in the opposite direction.

**Change.** When the resolved effect type is `odds_ratio` (or `trait_type != "quantitative"`),
refuse `use_ncontrol_as_total` instead of assuming. Separately, add an explicit `n_col` /
`neff_col` mapping so total N and effective N are never carried in the same field.

### C3 — An explicitly configured `pvalue.type` bypasses the negative-value guard, turning signed `log10(p)` into p = 1

`p_values.py:397-400` returns before the safety check:

```python
configured = policies.pvalue.type
if configured != "auto":
    evidence["decision_source"] = "config (pvalue.type)"
    return configured, evidence          # <- returns before the min < 0 refusal
```

The `minimum < 0.0` refusal at `:406-414` and the median sanity check run **only** on the
`auto` path. A user who declares `pvalue.type: neglog10` on a column that is actually
`log10(p)` (negative for every genome-wide-significant variant) reaches
`convert_negative_log10_to_p_value`, whose floor `pvalue.mlogp_min` defaults to `0.0`
(yaml:3093). Every negative value is floored to 0.0 — i.e. **p = 1.0** — with only a
`ctx.qc(..., warn=True)` line.

**Why it matters.** ~100 % of variants get p = 1. SE is then derived from that p at step 08 and
Z at step 09, so the entire effect-statistics block is fabricated, and the run completes with
status OK.

**Change.** Run the `minimum < 0` refusal and the median sanity check on **every** path. A
configured scale is a declaration to be validated against the data, not a reason to skip
validation — which is the module's own stated design principle ("declared, not guessed…
ambiguous or incompatible inputs are refused rather than guessed", README:95-97).

### C4 — `info.out_of_range: clip` (the default) turns any out-of-range value into a *perfect* imputation score

`imputation_quality.py:444-448`:

```python
if out_of_range == "clip":
    df = df.with_columns(pl.col(info_col).clip(clip_min, clip_max).alias(info_col))
```

`out_of_range_mask` is built against `tolerance_limit` (`clip_tolerance: 1.05`) but the clip
itself uses `clip_max` (1.0). Legitimately, this rescales `r2hat` overshoot of 1.02. But the
branch is unbounded: a mis-mapped column — INFO = 87, or a `−log10 p` column mapped to
`imp_info_col` — becomes `INFO = 1.0`, the best possible score, for every variant.

Harmonisation itself never filters on INFO (`info.low_quality_threshold` is used only for
counting, `imputation_quality.py:353-396`) — which is the right design — so the corrupted value
propagates into the exported VCF and the `filtering` / `qc_summary` modules downstream then
pass everything.

**Change.** Clip only within `[clip_min, tolerance_limit]`. Anything beyond the tolerance is
not rounding error: send it to `null`/`reject`, or fail on a fraction-based guard. The three
independently re-implemented "bounded quantity with rounding tolerance" behaviours
(`pvalue.tolerance_above_one` → reject, `eaf.clip_tolerance` → reject, `info.clip_tolerance` →
clip-everything) should share one helper; that alone removes this class of bug.

### C5 — Polars thread pool is never bounded: up to 1000+ runnable threads per run

`_derive_parallelism` (`service.py:1159`) computes `(workers, threads)`, but `threads` is
passed **only** to `annotate_and_liftover_vcf(threads=…)`, i.e. to bcftools. A repo-wide grep
for `POLARS_MAX_THREADS` / `set_num_threads` / `thread_pool_size` returns **zero hits**, so
every spawned worker builds a Polars pool of `os.cpu_count()`.

With `--threads 16` on a 64-core host: 16 workers × 64 Polars threads = **1024 runnable
threads**, while bcftools — the one consumer actually told what to use — gets `16 // 16 = 1`.
The unbounded consumer is uncapped and the bounded one is starved.

**Change.** Pass an initializer to the worker pool that sets
`os.environ["POLARS_MAX_THREADS"] = str(threads)` **before** `import polars` in the child
(spawn context), and let `total_cpu_budget` set `workers = budget // threads_per_chromosome`
rather than collapsing `threads` to 1. Related: `core/execution/runtime.py:171
auto_detect_ram_gb` and `:250 safe_thread_count` have **zero callers** — worker count never
scales with available RAM. Wire `safe_thread_count` into `_derive_parallelism` or delete both.

### C6 — Chromosome partitioning holds ~6 full-size copies of the dataset

`chromosome_partition.py:81-119`:

```python
source_for_partition = source_snapshot.join(assignments, on=..., how="inner")   # :81
... source_for_partition.partition_by(partition_col, as_dict=True, ...)          # :91
partitions = df.partition_by(chr_col, as_dict=True, maintain_order=True)         # :109
```

At this moment the process holds `df`, `source_snapshot` (every original column, every row),
the full-size join result, **and** two `as_dict=True` partition dicts that each materialise
every partition simultaneously. For a 50 M × 15 sumstats that is roughly six live full-size
copies — this is where 10–100 M-row datasets fail.

**Change.** Loop `for chrom in observed_chromosomes:` and write one chromosome at a time —
`df.filter(pl.col(chr_col) == chrom)` and `source_snapshot.join(ids_for_chrom, how="semi")` —
releasing each subset before the next. Drop the whole-frame `source_for_partition` join; the
row-count invariant it enforces is checkable per chromosome from the same counts. Peak RSS
falls from ~6× to ~2× the frame. Related and nearly free: `subset.write_csv` (`:118`) →
`write_parquet`, with `pl.read_parquet` at `service.py:1543`. The source snapshot beside it is
already Parquet; the working table pays a full text serialise, 3–5× the disk, a full re-parse
and a schema-inference scan for data whose schema is known exactly.

---

## Important

### I1 — Reconstructed BETA/SE on chromosome X assume autosomal variance

`effect_from_z.py:693-700` and `:749-756` build `2·EAF·(1−EAF)·Neff` with no chromosome
branch; `chromosome` appears only in error strings. X is in both `chromosome.allowed`
(yaml:1655) and `chromosome.allowed_after_split` (yaml:1487) — only MT and Y are dropped.
Repo-wide, nothing mentions hemizygosity or dosage compensation.

**Exposure is narrow but real:** this is Combination 4 of 4 in `effect_from_z.py:512`, reached
only when the study supplies **neither** BETA nor SE and does supply Z. The other three
combinations (`:348`, `:400`, `:462`) are exact and frequency-free.

**Change.** Either exclude X from Z-only reconstruction, or add an explicit X policy (male
fraction / `Neff_X`) and state the assumption. At minimum, emit a per-chromosome warning when
X reaches step 06 via Combination 4.

### I2 — `harmonise_effect_estimates` corrupts SE when the OR column is literally named `beta`

`effect_type.py:816` and `:832-837` write `.alias("beta")`. If the configured `beta_or_col` is
itself named `beta` — entirely possible, since effect type is auto-detected from the data
distribution, not from the column name — the OR is replaced in place by `log(OR)`. The
*subsequent* `with_columns` at `:870-878` then evaluates `pl.col(se_col) / pl.col(effect_col)`
with `effect_col` now holding `log(OR)`: SE is divided by the log-odds instead of the OR, and
is silently **nulled for every variant with OR < 1** (because the `when(effect_col > 0)` guard
now tests `log(OR) > 0`). Only `SOURCE_INPUT_ROW_COLUMN` is reserved
(`summary_statistics_io.py:1427-1431`). Same class of collision: `effect_from_z.py:406,469`
(`SE`, `BETA`), `z_score.py:137`, `sample_size.py:400,427,573`.

**Change.** Reserve these output names in the guard at `summary_statistics_io.py:1427`, or
allocate them through the existing `_unused_name` helper (`shared/allele_join.py:12-17`).

### I3 — Duplicate detection cannot see reverse-complement duplicates

Dedup runs at the dataset stage (`summary_statistics_io.py:1660`) on a `chr,pos,ea,oa` key
with lexicographic allele normalisation (`:1062-1089`) — so `A/G` and `G/A` collapse
correctly. But `A/G` and `T/C` at one coordinate — the same variant reported on opposite
strands, which happens in merged and meta-analysed files — form two distinct keys, survive
dedup, and are both flipped to `A/G` at step 04. The exported VCF then contains a true
duplicate and nothing checks for it.

**Change.** Add a post-orientation duplicate check at step 12/13, or re-key on the
strand-canonical allele pair, rejecting with `duplicate_variant`.

### I4 — `validation.max_reject_fraction` defaults to 1.0, so the rejection-fraction guard can never fire

yaml:4215 sets the default to `1.0`; only `warn_reject_fraction: 0.2` (yaml:4189) acts. The
YAML's own help text concedes "the test can never fire". The docs present this as a protective
gate (`processing-order.md:464`).

**Change.** Ship the registry's own recommended `0.95`, or state plainly in the docs that
harmonisation warns above 20 % removal and never fails at the shipped default.

### I5 — The `-log10 p` branch skips the range and missingness gates

`p_values.py:1117-1150`: the `raw` branch calls `_harmonise_raw_pvalues`, which honours
`pvalue.out_of_range` (`:729-743`) and rejects nulls as `pval_null`. The `neglog10` branch
calls only `convert_negative_log10_to_p_value` and applies neither. `n_invalid` is counted at
`:543` and stored as `mlogp_invalid_or_missing` (`:656`) but acted on nowhere. So
`pvalue.out_of_range: fail`/`reject` are silently inert for `-log10` studies, and null/`inf`
LP rows survive to step 13, where they are attributed to `final_missing_pval` instead of
`pval_null` — reject-reason attribution differs by input representation.

**Change.** Route both branches through one shared range/missingness gate applied *after*
conversion.

### I6 — Dataset step 3 rejects coordinates and alleles twice, and the on-screen ledger never reconciles

`summary_statistics_io.py:1531` runs `harmonise_coordinates_and_alleles`, which already rejects
null chromosome, null/sub-minimum position and disallowed chromosomes. The "CHR FIX" block at
`:1588-1610` then re-applies the *same* predicates and computes
`removed_coords = before_chr − df.height`, which is therefore ~always 0. That value is printed
as "Invalid coordinates" (`service.py:683-735`, `:2247`), so
`Input − missing − invalid_coords − non_standard − duplicates ≠ Ready for harmonisation` for
any file containing Y/MT or unparseable positions. The reject *file* is correct; the
human-facing ledger is not.

The second pass is not entirely dead — a position that parses as text but fails the numeric
cast at `:1563-1578` is caught only there. That is the documented order's real defect (see
[Documentation, §Order](#order-discrepancies)).

**Change.** Take the ledger counts from `RejectCollector.counts()` by reason rather than
recomputing frame deltas; keep the post-cast pass but scope it to post-cast failures only, and
document it as a distinct check.

### I7 — Reference panels are read whole, then projected — in every worker concurrently

`core/io/tables.py:15` and `shared/variant_columns.py:200,262`: `pl.read_csv(path, …)` /
`pl.read_parquet(path)` followed by `.select(required_columns)`. A per-chromosome 1000G/gnomAD
frequency table with 20+ population columns is fully materialised so five columns can be kept.

**Change.** Push the projection into the reader — `pl.read_csv(..., columns=required_columns)`
and `pl.read_parquet(path, columns=required_columns)`. The "absent column" check moves ahead of
the read, which `read_delimited_table` already effectively does via `resolve_delimiter`.

### I8 — Two failure paths destroy the evidence needed to diagnose them

**(a)** `service.py:1118 _reject_counts_from_file` catches bare `Exception` and returns `{}`. A
corrupt or unreadable reject file makes the `input` column of
`qc_summary/{sample}_reject_reasons.tsv` read `0` for every reason — indistinguishable from
"nothing was rejected", in a table whose own docstring says a zero is *positive evidence that
the check ran*. Everywhere else in `harmonise_chromosomes`, missing reject provenance is a
hard `PipelineError`. Use `pl.scan_csv(...).group_by(...).len().collect()` and let the
exception surface as a `PipelineError` naming the path.

**(b)** `service.py:1958 except BaseException` records only `"%s: %s" % (type(exc).__name__, exc)`.
The comment assumes `logger.step()` already logged the traceback, but the code *between* steps
(the orientation-decisions block `:1670`, the reconciliation `:1847`) is inside no step at all,
so a `KeyError: 'chr_col'` reaches the user as exactly that string — no file, no line, no
context beyond the chromosome name. It also catches `MemoryError` and returns
`status="failed"`, so an OOM is retried up to `max_retry_rounds` under identical memory
pressure. Add `qc_dict["traceback"] = traceback.format_exc()` and re-raise `MemoryError`.

### I9 — Three divergent missing-value vocabularies

`core/values.py:9-11` (`"", NA, N/A, NAN, NONE, NULL, -, ., NOT_APPLICABLE`) governs scalar
config parsing; `input.null_values` (`['NA','NAN','','.']`) governs the CSV parse
(`summary_statistics_io.py:124,163`); `input.chromosome_null_values`
(`['','NA','NaN','.','null']`) is a third list with different casing. Harmless for numeric
columns (non-strict cast → null) but **not** for identifiers: `variant_identifiers.py:41-56`
derives its sentinels from `input.null_values`, so a literal `-` is treated as a present rsID
and exported as one.

**Change.** One shared `missing_tokens()` used by the reader, the identifier step and
`optional_text`.

### I10 — Indels are never normalised before the reference join

`strand.py:484-489` correctly restricts reverse-complement to 1 bp alleles, and `I`/`D`/`-`
codings are rejected upstream by `allele.pattern`. But the join is CHROM/POS only
(`strand.py:466-470`) followed by exact string equality, and no left-alignment or
padding-base trimming exists anywhere in the module — `bcftools norm` runs only downstream
(`vcf_processing.py:509-525`). A right-aligned or differently-anchored indel silently becomes
`reference_unmatched` and is removed under the default `strand.unmatched_action: reject`.

**Change.** Left-align and trim study and reference indels against the FASTA before the join,
and split the `reference_unmatched` counter into SNV vs indel so the loss is visible.

### I11 — Whole-dataset reject concatenation, then a third pure-Python pass

`rejects.py:707-770 concat_reject_files` reads every per-chromosome reject file into `frames`
and combines with `pl.concat(frames, how="diagonal")` — all 23 frames plus the result live
simultaneously, with *all original columns*. Then `_validate_reject_output` (`:668`) re-reads
the written (possibly gzipped) file with `csv.reader` and `sum(1 for _row in reader)` — a
third full pass in pure Python, despite the docstring claiming it works "without loading it
back into memory". On a dataset where the read gate rejects 30 % of 80 M rows this is ~24 M
rows held twice.

**Change.** `pl.scan_csv(paths).sink_csv(temporary)`; take the row count from
`pl.scan_csv(temporary).select(pl.len()).collect()` and the header from the first line only.

---

## Optional

- **O1 — `_inspect_text_file` decompresses the whole input in a Python `for` loop**
  (`summary_statistics_io.py:1753`) purely to count lines and test whether the last line ends
  with `"."`. Replace with a chunked binary count (`while chunk := f.read(1<<20): n += chunk.count(b"\n")`,
  keeping the tail for the truncation test) — ~15× faster, 1–4 min of single-threaded CPU
  saved per dataset before any work begins.
- **O2 — a full extra copy of every raw chromosome VCF.** `vcf_processing.py:460`
  `shutil.copy2(input_vcf, original_vcf)`; `_ORIGINAL` is used only for
  `counts["ORIGINAL"] = get_vcf_variant_count(...)` (`:494`) — a number `run_gwas2vcf` already
  verified against `expected_variants` (`gwas2vcf_runner.py:126`) — and is deleted again at
  `:1332`. Take the count from `input_vcf` directly.
- **O3 — a QC audit file is rewritten in place by a function that only claims to read it.**
  `service.py:3594-3609` does `pd.read_csv(gwas2vcf_summary_path)` … `to_csv(same path)`. A
  partial write truncates the audit trail, and the `except` reports only "cannot prepare the
  pre-VCF QC counts". Strip the repeated header rows in `cleanup._concatenate` where the file
  is built. (Also the last pandas dependency in the hot path, alongside `_save_qc_results` /
  `_metric_total`.)
- **O4 — out-of-range INFO can be rejected under the wrong reason.**
  `imputation_quality.py:459-474` blanks out-of-range values under `info.out_of_range: null`;
  `:496-498` folds them into `missing_now` and `:510-524` rejects them as `info_missing`, even
  though `rejects.py:184-187` declares `info_out_of_range` and `info_missing` as distinct,
  mutually exclusive reasons.
- **O5 — a literal `p = 0` is counted as a "positive p-value below Float64".**
  `p_values.py:224-236` — the underflow append is outside the `if/elif` chain, so a token of
  `"0"` lands in both `reported_zero_indices` and `underflow_indices`, propagating to
  `initial_positive_pvalues_below_float64` and `pvalue_float_underflow`. The scientific path is
  still correct (`exact_underflow` at `:700-704` additionally requires a finite `LOG_P`); only
  the QC metrics are wrong — and they contradict the module docstring's own claim. Make the
  flag `elif`-exclusive.
- **O6 — `Neff = 4/(1/Ncase + 1/Ncontrol)` has no zero guard** (`sample_size.py:69-73`). With
  `Ncase = 0`, Polars gives `1/0 → inf` and `Neff → 0`, which with `sample_size.min_value: None`
  reaches step 06 as a legitimate-looking zero. Gate on `cases > 0 & controls > 0`.
- **O7 — dataset-level row accounting is not a hard assertion.** `service.py:3134-3148` sums
  over `completed` chromosomes only and reports `balanced` as a boolean that is tautologically
  true (the per-chromosome `reconcile()` at `:1854` already raised otherwise). Rows in a failed
  chromosome, and every row removed at the dataset stage, are outside the identity — so
  "input = output + rejects" is never verified end to end.
- **O8 — `position.min_value: 0` permits position 0**, which is not a valid 1-based coordinate,
  while `rejects.py:76-77` describes the check as "not a positive number". Set the default to 1.
  Related: `shared/variant_columns.py:160-161` casts position `Float64 → Int64`, so `"100.7"`
  silently truncates to 100 and `"1e5"` becomes 100000; only the *split* path
  (`coordinates.py:394-457`) applies `position.extraction`.
- **O9 — `SE = BETA / Z` is not absolute-valued** (`effect_from_z.py:305`, `:403-406`).
  `usable_for_se` only checks `Z != 0` and finiteness, so a sign-inconsistent study Z yields a
  negative SE, caught only at `effect_validation.py:300-307`. `standard_error.py:647` correctly
  uses `np.abs`.
- **O10 — non-positive OR keeps swapped alleles but not the flip.** `strand.py:792`
  `pl.when(swap & (effect > 0.0)).then(1.0/effect).otherwise(effect)` — for `swap=True` with
  `effect ≤ 0`, alleles are rewritten (`:786-788`) while the effect is returned unchanged, with
  no reject reason. Unreachable through `service.py` today (step 03 always normalises to
  log-odds first), but `harmonise_strand_orientation` is in `__all__`. Use `.otherwise(None)`.
- **O11 — reverse-complement matching is applied to multi-base alleles in build detection.**
  `genome_build.py:268-271` builds `reverse` without the SNV guard defined at `:272-277`; the
  `informative` counts are guarded (`:292-301`) but `matched` (`:285-288`) is not, so an indel
  "reverse-complement match" can inflate the build-selection statistic.
- **O12 — dead code and inert config.** `service.py:193` exports
  `"extract_chromosome_from_filename"` in `__all__`, which **does not exist anywhere in the
  repo** — `from …service import *` raises `AttributeError`.
  `summary_statistics_io.py:1805 check_file_truncation` has zero callers.
  `sample_size.cases_only` is declared (yaml:3692), listed in the module's policy list
  (`sample_size.py:46`) and named in a hard-coded error *string* (`:597`), but never read —
  though its enum is single-member `[fail]`, so no user-reachable value diverges from the
  hard-coded behaviour. `execution.threads_per_chromosome` (default 5) is unreachable through
  the CLI: `cli.py:493` always sets `total_cpu_budget = config.execution.threads`, and
  `_derive_parallelism` then overwrites `threads` with `budget // workers`.
- **O13 — multi-cohort INFO is summarised by an unweighted row-wise median**
  (`summary_statistics_io.py:420-438`) where a sample-size-weighted mean is standard, is
  reported via `print()` rather than the logger, and runs *after* dedup (`:1709` vs `:1660`) —
  so INFO-based duplicate tie-breaking (`:1228-1239`) is inert for exactly those studies.
- **O14 — step-label drift in reject provenance.** `imputation_quality.py:1` says "Step 12" but
  `STEP_LABEL = "11 info_harmonisation"`; `variant_identifiers.py:1` says "Step 13" but
  `STEP_LABEL = "12 snp_column"`.

---

## Consolidation opportunities

Listed only where the copies have already diverged in behaviour or where the shared version
fixes a defect above.

| Duplicated logic | Locations | Proposal |
|---|---|---|
| **Partition a delimited table to per-chromosome Parquet** | `concordance/service.py:208,355` (lazy, `sink_parquet` — correct); `external_reference_staging.py:139,206` (batched `ParquetWriter` — correct); `chromosome_partition.py:14` (fully eager — the one that OOMs) | One `partition_table_to_parquet(source, destination_for, *, partition_expression, columns, partitions, separator, null_values, infer_schema_length, compression, atomic_suffix, error_type) -> dict[str, int]` in `core/io/tables.py`. Highest-value item here: it fixes **C6** by construction and gives the reject/snapshot writers the atomic-rename discipline `external_reference_staging.py:266-290` already has. |
| **Chromosome sort key** | `service.py:738` (`X→100, Y→101, XY→102, MT→103`) vs `resource_preflight.py:49` (`(1, text)` for all non-digits → `MT, X, XY, Y`) | Same name, **different orderings** — preflight failure reports list chromosomes differently from every other output. Move the service version to `shared/runtime.py`. |
| **"Is this orientation action swapped?"** re-derived by string matching | `allele_frequency.py:511-513, 528-533, 817, 1277-1279` — while the authoritative boolean `__strand_swapped` (`strand.py:513`) is *dropped* at `strand.py:819-826` | Export `__strand_swapped`, or a shared `strand_action_is_swapped()`. |
| **Forward / swap / RC / RC-swap candidate enumeration + SNV guard** | `strand.py:484-507` vs `genome_build.py:263-278` — two independent implementations, and the second is missing the SNV guard (**O11**) | One helper. |
| **Valid-frequency mask** | `shared/statistics.py:8-20` plus hand-rolled copies at `strand.py:252-263`, `allele_frequency.py:323-330, 866-871, 1044-1049` — which also differ in whether they cast to Float64 first, so a Utf8-inferred frequency column raises `ComputeError` in some branches and not others | Use the shared one everywhere. |
| **Bounded quantity with rounding tolerance** | `pvalue.tolerance_above_one` (reject), `eaf.clip_tolerance` (reject), `info.clip_tolerance` (clip-everything) | One helper; removes **C4**. |
| **Degenerate-EAF (0/1) rejection** | `allele_frequency.py:1572-1601` and again `effect_from_z.py:612-660` | One helper. |
| **Numeric-with-separators parsing** | `sample_size.sample_count_expression:61-66` and `summary_statistics_io._numeric_ranking_expression:975-980` — the same `[,\s_]`-stripping cast, written twice | One helper. |
| **Safe output-name check** | `service.py:746 _safe_glob_name` vs `cleanup.py:55 _safe_name` (one also rejects `.`/`..`) | `shared/runtime.safe_output_name(value, what, error_type)`. |
| **`_emit`** | Redefined four times: `cleanup.py:52`, `summary_statistics_io.py:71`, `input_validation.py:290`, and `vcf_processing.py:81` with an **incompatible** `(message, logger, level)` signature | Align `vcf_processing._emit` with `shared.runtime.emit_message` and drop the rest. |

### Function size

Measured: `harmonise_chromosomes` 591, `process_one_chromosome` 562,
`_prepare_dataset_for_chromosome_processing` 559, `run_harmonisation_pipeline` 437,
`_run_post_merge_stages` 418.

- **Leave alone.** `process_one_chromosome` and `_run_post_merge_stages` are genuinely linear —
  16 and 5 numbered steps, no branching. Splitting them would only add parameter passing.
- **Worth doing.** The `_..._block` display functions plus `_announce`/`_count`
  (`service.py:273-737`) are 449 lines of pure `dict → str` rendering, identical in shape to
  `qc_reporting.harmonisation_qc_takeaway_lines`. Moving them into `qc_reporting.py` is
  mechanical and behaviour-preserving. `_prepare_dataset_for_chromosome_processing` is ~40 %
  screen rendering (`:2308-2360`, `:2504-2537`, `:2545-2595`) and shrinks accordingly.
  `harmonise_chromosomes` mixes the retry loop with ~150 lines of reject-file reconciliation
  (`:3010-3140`) that has nothing to do with fan-out — extract
  `_finalise_dataset_rejects(...) -> dict`.
- **Nothing else** in the restructure bucket earns its cost.

---

## Documentation review — harmonisation only

The README's harmonisation section (`README:67-129`) is short by design and delegates to
`docs/wiki/harmonisation/processing-order.md`. That is the right structure, and
**processing-order.md's stage and step *sequence* (`:69-100`) matches the code exactly.** The
problems are inside step descriptions, not in the order.

### Critical

| Finding | Evidence | Corrected statement |
|---|---|---|
| **Every supplied variant ID has `:` rewritten to `_`.** A study with IDs `1:12345:A:G` exits as `1_12345_A_G` and cannot be joined back by ID. Docs say only "Preserve and normalize a supplied study ID" (`processing-order.md:471`) | `variant_identifiers.py:197-207` | "Step 12 rewrites every colon in the identifier to an underscore so the ID survives the VCF round trip; the count is reported as `identifiers_rewritten`." |
| **`field_lifecycle.snp.note` in the YAML says IDs are "Constructed as chr:pos:ea:oa"** — and this string is surfaced to users in the field-completeness TSV and manifest | `_identifier_expr` builds `chr_pos_ea_oa` with underscores (`variant_identifiers.py:68-79`) | Change the note to `chr_pos_ea_oa`. Also `field_lifecycle.snp.recovered` says "chromosome step 13"; identifiers are built at **step 12** (`service.py:1830`). |
| **`eaf.degenerate: reject` (default) removes every variant whose EAF is exactly 0 or 1 — for all variants**, not only Z-reconstruction candidates. Docs mention the (0,1) requirement only for Z-only recovery (`processing-order.md:431-433`) | `allele_frequency.py:1100, 1571-1576` | "Step 4 rejects EAF of exactly 0 or 1 for every variant under `eaf.degenerate` (default `reject`), reason `eaf_degenerate`." |
| **The reject-reason table's `step` column contradicts the run.** `rejects.REASON_STEPS` is a legacy numbering that disagrees with both the runtime step numbers and the `reject_step` label written into the row-level file: `eaf_*`/`palindromic_*`/`af_discordant`/`reference_*` → table says **03**, runtime is **04**; `neff_invalid`/`sample_size_invalid` → table **04**, runtime **05**; `effect_non_positive_or` → table **05**, runtime **03** | `rejects.py:229-241` vs `allele_frequency.py:89`, `sample_size.py:41`, `effect_type.py:61` | Renumber `REASON_STEPS` to the runtime steps. (`outputs-and-qc.md:101-102` does not mention the column at all.) |
| **Z-only reconstruction produces a standardized effect, not a log-odds ratio, for case-control studies** — despite step 03 having just placed supplied effects on the log-odds scale. Docs say only "a standardized effect estimate" (`processing-order.md:588`) and never name the method | `effect_from_z.method` default `metal_large_n` (alt `zhu_2016`); `phenotype_standard_deviation` default `null` (σ_y = 1) | Name both config keys and state the scale explicitly. This is the single most consequential thing a user can misread. |

### Important

- **`pvalue.out_of_range` defaults to `clip`** (floor `pvalue.clip_low: 1e-300`), so
  `pval_out_of_range` never fires under packaged settings — and any SE derived from a clipped p
  is then **rejected** at step 10 (`validation.se_from_clipped_pval: reject`,
  `effect_validation.py:471-478`). None of this is documented.
- **`sample_size.missing_action` defaults to `remove`** with `max_missing_fraction: 0.01`.
  `processing-order.md:272-276` mentions only mean/median filling and "too much missing … stops
  the dataset". Variants with missing N are dropped by default.
- **`final_check.require` is `[chr, pos, eaf, beta, se, zscore, pval]`** — `snp`, `ea`, `oa`,
  `info` and `n` are deliberately **not** gated at step 13. Docs say only "every required
  column" (`processing-order.md:476`).
- **`--zero-p-se-action` (`fail`/`approximate`/`reject`) appears in no document** (`cli.py:472`).
- **`processing-order.md:464` presents `max_reject_fraction` as a protective guard**; at the
  shipped default it can never fire (see [I4](#i4)).
- **Three produced outputs are never named**:
  `qc_summary/{id}_{build}_population_frequency_qc.json`,
  `qc_summary/{id}_field_completeness.tsv`, and
  `{id}_notlifted_{target_build}_merged.vcf.gz`.
- **The FinnGen resource key is spelled `fingen`**, not `finngen`, in both the YAML
  (`available_sources`, yaml:38,44) and the README (`README:414`). The documentation is
  *accurate*, but the spelling will cost users a support round-trip; rename the key (with an
  alias) or call it out explicitly.
- **`processing-order.md:7-9` claims "the same stage and step names appear in the log files".**
  The step *numbers* match; the titles are paraphrases (step 5 logs "Genome build" vs doc
  "Determine the genome build"; step 7 "P-value scale" vs "Harmonise P values"; step 11
  "Imputation quality (INFO)" vs "Obtain INFO"; step 12 "SNP identifier" vs "Complete variant
  identifiers"). Say the numbers match and the titles are paraphrased.

### Order discrepancies

One, inside dataset step 3 (`processing-order.md:232-249` claims a strict "in order" list):

| Doc order | Code order (`summary_statistics_io.py:1462-1702`) |
|---|---|
| 1 missingness gate | 1 missingness gate (`:1462-1516`) |
| 2 split combined coordinate | 2 `harmonise_coordinates_and_alleles` — split, normalise, **and reject** (`:1528-1537`) |
| 3 normalise chr/pos/alleles | 3 numeric casting (`:1563-1578`) |
| 4 reject invalid coordinates/alleles | 4 **second** rejection pass: `invalid_chromosome`, `invalid_position` (`:1581-1610`) |
| 5 numeric conversion | 5 **second** allele pass: `null_allele`, `non_standard_allele` (`:1613-1643`) |
| 6 duplicates / 7 internal INFO | 6 duplicates (`:1655`) / 7 internal INFO (`:1700`) |

The user-visible consequence: a position that parses as text but fails the numeric cast is
rejected as `invalid_position` in the *post-cast* pass, which the documented single-pass
ordering cannot explain. Fix the code first ([I6](#i6)), then document the surviving check.

### The main gap for a user tracing their data end to end

**There is no consolidated table of default dispositions.** Nowhere can a user see in one place
that the shipped policy is `reject` for `eaf.out_of_range`, `eaf.degenerate`,
`validation.se_invalid`, `validation.beta_invalid`, `validation.z_invalid`,
`validation.se_from_clipped_pval`, `effect_from_z.low_effective_variance_action` and
`sample_size.missing_action`; `keep` for `validation.beta_zero` and `info.on_missing`; `clip`
for `pvalue.out_of_range` and `info.out_of_range`; and `warn` for
`strand.af_discordance_action`, `validation.beta_se_z_concordance` and
`validation.z_pval_concordance`. `processing-order.md:460-469` says "follow the resolved
policy" without naming a single default. **This is the highest-value documentation change in
this review** — one table, roughly 20 rows, generated from the YAML.

Secondary gaps:

- **Recomputed vs preserved is never stated per field.** The `field_lifecycle` block in the YAML
  holds exactly this information (and is partly wrong — see Critical above) but is surfaced in
  no document.
- **The 48 reject reasons carry plain-English descriptions** (`rejects.py:67-215`) that are
  written into the `description` column, yet no document lists or links them.
- **Harmonisation does not filter on INFO.** `info.low_quality_threshold` is used only for
  counting (`imputation_quality.py:353-396`); filtering happens in the `filtering` module. This
  is the correct design and should be stated, because users reasonably assume the opposite.
- **`PARTIAL` dataset status is not explained.** `_combined_dataset_status` (`service.py:1050`)
  returns it from either chromosome failure or merge failure; the docs name the label but not
  that the merged VCF is then missing whole chromosomes.

### Verified as accurate

The "six basic SE, BETA and Z checks" (`effect_validation.py:435`);
`Neff = 4/(1/Ncase+1/Ncontrol)`; the PLINK `23/24/25/26/XY/PAR1/PAR2/M` rename map; the three
promised primary VCFs (`service.py:1021-1047`); all seven `resource_layout` templates and the
ALFA/1000G/wgs_ukb/panukb/finngen source lists (`README:389-417`); the 27-column sample-sheet
schema with `extra="forbid"` (`sample_sheet.py:148-176`); and every one of ~25 spot-checked
config keys and defaults (`build.min_reference_match_fraction` 0.5,
`strand.consensus_threshold` 0.99, `strand.min_informative_variants` 1000,
`pvalue.zero_missing_se` `fail`, `pvalue.se_tail` 2, `validation.se_division_floor` 1e-12,
`effect_from_z.minimum_effective_variance` 1.0, `vcf.liftover_swap` `exclude`,
`vcf.on_merge_failure` `fail`, `duplicates.conflicting_action` `remove_all`, …). No documented
key was missing or renamed.

---

## Verified correct

Noted only because a reviewer might otherwise suspect them.

- **Combined complement + swap** — the classic sign trap. `reverse_complement_swapped` sets
  `__strand_swapped=True` and `reverse_complement` sets it `False` (`strand.py:493-506`).
  Correct.
- **Exactly one flip, never two.** β, Z and EAF all use the same `swap` mask
  (`strand.py:784, 794, 800, 803`); SE and p are correctly untouched. The MAF path defers the
  EAF flip to `allele_frequency.py:1276-1285` precisely because `strand.py:801` skipped it —
  so no double flip on any path.
- **Complementing** uses a single-pass `replace_many` plus `str.reverse()`
  (`shared/variant_columns.py:83-90`), so no A→T→A chaining; non-ACGT is nulled.
- **Palindromic policy is stricter than common practice.** Informative band [0.40, 0.60]
  (yaml:2288, 2316) with a 0.20 minimum error margin (yaml:2361) — versus the usual 0.42
  threshold. Unresolvable palindromes are *always* dropped with an explicit reason
  (`strand.py:743-777`; the action enum is `reject|fail` only), and the candidate frequency
  compared to the panel is correctly re-oriented to the reference ALT (`strand.py:546-552`).
- **Multi-allelic sites** leave more than one candidate and are rejected as
  `reference_ambiguous` (`strand.py:653-657, 778-782`) — not silently treated as duplicates.
- **Derived-statistic formulas are algebraically correct and internally consistent.**
  `BETA = σ_y·Z/√(2p(1−p)·S)` and `SE = σ_y/√(2p(1−p)·S)` with `S = Neff` or `Neff + Z²`
  (`effect_from_z.py:749-758, 835-849`) reproduce `Z = BETA/SE` exactly and match METAL / Zhu
  et al.; `SE[logOR] = SE[OR]/OR` (`effect_type.py:370`); `Z = β/SE` (`z_score.py:135`).
- **Extreme p-values are handled in log space** via `log_ndtr` (`shared/statistics.py:23-32`)
  and `ndtri_exp` (`core/statistics.py:41-55`), so p < 1e-308 does not underflow.
- **Numerical guards are present where they matter**: EAF endpoints rejected before Z
  reconstruction (`effect_from_z.py:615-617`), denominator positivity re-checked per row
  (`:760-778`), `se_division_floor = 1e-12`, and finiteness tested *before* positivity so `inf`
  cannot pass a `> 0` test (`effect_validation.py:293-307`).
- **No OR-confidence-interval path exists** anywhere in the module — SE must be supplied, or
  derived from Z or from p. This is a deliberate gap, not a bug, but it is worth documenting so
  it does not read as an oversight.

---

## Suggested order of work

1. **C3, C4** — one-line-ish policy/guard fixes with immediate scientific payoff.
2. **C1, C2** — the two missing cross-checks. Both use evidence already in memory.
3. **C6 + the `partition_table_to_parquet` consolidation** — one change that fixes the memory
   ceiling and removes the third copy of the partitioning logic.
4. **C5** — worker/thread initializer; independent of everything else.
5. **The default-disposition table in `processing-order.md`**, plus the four Critical doc
   corrections. Cheap, and it is what users actually need.
6. Important items as capacity allows; Optional items opportunistically.
