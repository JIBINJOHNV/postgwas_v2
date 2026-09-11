# PostGWAS `ld_clump` — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/ld_clumping/` (`ld_prune_standard.py` 2822, `service.py` 632, `ld_prune_region.py` 285, `common.py` 151, `cli.py` 132), `config/models/modules/ld_clumping.py`, `config/defaults/modules/ld_clumping.yaml`, `docs/wiki/modules/ld-clumping.md`, and the reference-layout claims in `examples/fine_mapping/README.md`, `README.md` and `docs/wiki/modules/flames.md`.
**Method:** read-only on the tree — no file was modified. The FUMA algorithm was traced stage by stage against the implementation. The performance claims in **I1** and **I2** are measured, not estimated: a faithful reproduction of `_inventory_lookup`'s validation loop was run against its vectorised equivalent, and the same loop was run under `ThreadPoolExecutor` to measure GIL scaling. Line numbers refer to the working tree as of 2026-08-26.

**Summary — and this is the headline: I did not find a scientific defect in the clumping algorithm.** The three FUMA stages are implemented correctly, including a subtle invariant that keeps the locus lead consistent with the locus's reported P (see *Verified correct*), and the reference-manifest validation in `service.py` is the strongest input contract anywhere in this codebase. The findings below are about **how the work is executed, not what it computes**: the chromosome parallelism cannot accelerate the algorithm's core, one validation routine is an order of magnitude slower than it needs to be, and — the one user-facing defect — the shipped worked example documents a reference layout that this module rejected two format versions ago, which under the shipped defaults yields a *successful* run reporting zero genomic risk loci.

There are no Critical findings.

---

## IMPORTANT

### I1 — Chromosome parallelism is thread-based, and the algorithm's core holds the GIL

**Where:** `ld_prune_standard.py:2539` (`ThreadPoolExecutor`), with the cost-aware scheduler and memory arithmetic at `:2466-2514`.

The module schedules chromosomes carefully — ordering by significant-variant count so dense chromosomes start first, and sizing workers from `gwas.estimated_size() × largest_share` against the memory budget — and then runs them on threads. But the dominant per-chromosome work is Python-level row iteration, which cannot run concurrently:

| routine | shape |
|---|---|
| `_inventory_lookup` (`:679`) | `for row in inventory.iter_rows(named=True)` over every LD endpoint |
| `find_ind_sig_snps` (`:1138`) | `for top_snp in ranked` — one greedy round per significant variant |
| `find_lead_snps` (`:1244`) | `for top_row in ind_sig_heads.iter_rows(named=True)` |
| `define_genomic_risk_loci` (`:1475`) | `for row in initial.iter_rows(named=True)` — the interval merge |

Measured on the `_inventory_lookup` loop (500,000 rows, four repetitions):

```
per-row loop ×4 tasks, 1 thread :  8.29 s
per-row loop ×4 tasks, 4 threads:  9.00 s
```

**No speedup — marginally worse.** What the threads *do* parallelise is real but partial: the `tabix` subprocesses and the polars calls inside each round release the GIL. So the scheduler is optimising a parallelism that applies to the I/O and the vectorised fragments, while the greedy clumping rounds and the inventory validation serialise.

**The tradeoff is not simply "use processes".** The current design deliberately hands the *whole* genome-wide frame to every worker and subsets inside (`:2544-2548`), with a comment explaining that pre-slicing in the submitting thread would materialise all 22 subsets before any worker starts. That is free under threads and catastrophic under processes, where each submission would pickle the entire study. Escaping the GIL therefore means either vectorising the loops first (I2 is the cheapest instance) or pre-partitioning per chromosome and paying one serialisation cost per worker. Worth deciding deliberately; right now the parallelism reads as broader than it is.

### I2 — `_inventory_lookup` validates row by row in Python; the vectorised form is 11× faster

**Where:** `ld_prune_standard.py:677-741`.

```python
for row in inventory.iter_rows(named=True):
    expected = canonical_variant_id(row["chromosome"], row["position"],
                                    row["allele_1"], row["allele_2"])
    valid = (... and row["canonical_id"] == expected)
```

Measured on 500,000 inventory rows, comparing this loop against the same five checks written as one polars filter:

| implementation | time | result |
|---|---|---|
| per-row Python loop | **2.19 s** | identical |
| vectorised polars filter | **0.20 s** | identical |

It is called **twice per chromosome** (`:773` on the index inventory, `:907` on the endpoint inventory), and the second call is the expensive one: `endpoint_positions` is the union of `pos_a` and `pos_b` across every retrieved LD pair, so its size scales with *verified indexes × partners* — the largest frame the module handles. The `break` after five invalid rows (`:699-700`) only helps when the reference is already broken; on a healthy reference the loop always scans everything.

The two duplicate-detection checks that follow (`:711-738`) are already vectorised polars aggregations. Only the row loop needs converting, and `canonical_variant_id`'s logic has a direct expression form — `add_canonical_ids` already writes it (`:161-181`).

### I3 — `find_lead_snps` does exactly what `find_ind_sig_snps`'s own comment forbids

**Where:** `ld_prune_standard.py:1258` versus the comment at `:1132-1133`.

`find_ind_sig_snps` explicitly avoids per-index rescans:

```python
# Walk the ranked candidates once and track depleted IDs in a set. Filtering
# the frame per index SNP rescans every remaining row and is quadratic.
```

`find_lead_snps` then does the rescan, once per lead:

```python
unassigned = ind_sig_heads.filter(~pl.col("uniq_id").is_in(list(assigned)))
members = unassigned.join(partners, on="uniq_id", how="inner")
```

Every iteration rebuilds `list(assigned)` from a growing set and re-filters the whole frame — quadratic in the number of independent significant SNPs, which for a well-powered GWAS is in the thousands. The fix is the pattern already written twenty lines up: keep `assigned` as a set, skip assigned heads on iteration (it already does, `:1246-1247`), and intersect `partners`' ids against the unassigned set in Python instead of re-filtering the frame.

### I4 — Every tabix payload is held in memory three times

**Where:** `ld_prune_standard.py:545-550`, then `:589` / `:629`.

```python
result = subprocess.run([tabix_bin, "--regions", ...], capture_output=True, text=True, check=True)
return result.stdout                      # 1. whole payload as a Python str
...
frame = pl.read_csv(StringIO(payload), ...)   # 2. StringIO copy   3. polars frame
```

For a dense chromosome the forward and reverse LD payloads are roughly *verified indexes × partners × 60 bytes*; at 5,000 verified indexes with 500 partners each that is ~150 MB per direction, and both directions plus the endpoint-inventory payload are live at once inside `prepare_chromosome_reference`. Because these run inside concurrent workers, the peak multiplies by the worker count.

Writing tabix's stdout to a temporary file (`stdout=handle`) and giving polars the path removes two of the three copies and lets polars stream. The module already uses `TemporaryDirectory` two lines above for the regions file, so the machinery is in hand.

### I5 — `missing_chromosome_action: warning` turns a wrong reference layout into a successful run with zero loci

**Where:** `ld_clumping.yaml:23`; `_missing_chromosome_references` at `:1944-2011`; summary at `:2270-2281`.

Under the shipped `format_version: 2` / `orientation: upper_triangle_dual_index` defaults, each chromosome requires **six** files: `{pop}_chr{n}.ld.gz`, `{pop}_chr{n}.reverse.ld.gz`, `{pop}_chr{n}.variants.tsv.gz`, and a `.tbi` for each. `_missing_chromosome_references` checks all six and, with the default action `warning`, records the chromosome as skipped and continues.

To the module's credit the reporting is honest and prominent — the summary heading changes to "Standard LD clumping finished with missing-reference warnings" and prints the skipped chromosomes with the number of significant variants dropped, `_write_reference_exclusions` persists them, and `service.py:584-596` propagates `status: "partial_reference"` into the pipeline context. What it does not do is fail. A reference directory in the wrong layout produces exit status 0, a report of zero genomic risk loci, and a warning line — and, as D1 below shows, the shipped worked example documents exactly that wrong layout.

---

## Documentation

### D1 — The worked example documents a reference layout this module no longer accepts

`docs/wiki/modules/ld-clumping.md` is accurate and current: it describes the forward table, the reverse sidecar, the variant inventory, the manifest, and the batched query strategy (`:48-52, :68, :169-173, :221-222, :296`).

The examples around it are not.

- `examples/fine_mapping/README.md:47-48` still describes the reference as "`reference/pairwise_ld/EUR_chr1.ld.gz` through the analysed chromosomes, each with a tabix index and seven columns", followed by the format-version-1 caveat that "every potential index SNP must therefore occur as variant A … Otherwise that SNP is explicitly logged as self-only". That caveat is precisely what the `reverse_file_pattern` sidecar was introduced to remove.
- `README.md:934` and `docs/wiki/modules/flames.md:109` point at the same `reference/pairwise_ld` layout.

Following that example gives every chromosome a missing `variants.tsv.gz` and a missing reverse sidecar, so `_missing_chromosome_references` skips all of them and — per I5 — the run *succeeds* with zero loci. There is also no `ld_reference.yaml`; without it `service.py:229-293` cannot validate format version, build, population, orientation, file patterns, column lists, `window_kb`, `minimum_r2` or `minimum_maf`, which is the module's best safety net.

Update the three example documents to the version-2 six-file-plus-manifest layout, and point them at the preparation tool that `ld-clumping.md:221-222` already names.

### D2 — Two locus definitions ship enabled by default with no guidance

`methods: [region, standard]` runs both. They answer different questions — `region` reports one lead per pre-computed LDetect block, `standard` performs FUMA clumping — and write two non-comparable sets of loci into the same output directory. The docs give the join key; they do not say which definition to carry into fine-mapping or reporting.

### D3 — Worth one line: `allele_mismatch` usually means the LD panel's strand, not the study's data

`canonical_variant_id` sorts alleles, so A/G and G/A are one variant, and palindromic SNPs match on either strand because their sorted pair is identical. A non-palindromic variant whose LD panel is on the opposite strand (study A/G, panel T/C) is reported as `allele_mismatch` and excluded. That is the conservative, correct behaviour and it is logged with `observed_reference_ids` — but a user seeing thousands of such exclusions should be told to check the panel's strand before they check their study.

---

## OPTIONAL

**O1 — `minimum_reference_maf: 0.00001` is not a MAF filter.** At 1e-5 it excludes monomorphic reference variants and nothing else. The name and the manifest check that guards it (`service.py:289-293`) both imply a real frequency floor. Either raise it to something scientifically meaningful or say in the YAML comment that its purpose is to reject monomorphic entries.

**O2 — the standard summary's `P_value` underflows to exactly 0.0 for any locus stronger than ~1e-308**, which is where lead SNPs live. `LP` is emitted alongside it, and the code is scrupulous about ranking on LP everywhere — but `ld_prune_region.py:194-196` carries an explicit "P is reported for readability only" comment and the standard summary does not.

**O3 — `pl.lit(False)` as a sort key** (`:1111`, `:1229`) is a placeholder that exists only so the two-element `descending=[False, True]` lines up in both branches. One comment saves the next reader a double-take.

**O4 — `_print_clumping_summary` prints unconditionally** (`:2229` onwards), like `sumstat_filter` and unlike `qc_summary`, which gates the equivalent on `configuration.logging.show_screen`.

**O5 — `candidate_p` in `define_genomic_risk_loci` is materialised as a Python dict** over every GWAS-tagged candidate (`:1317-1323`) purely to compute the three `summary_pvalue_thresholds` counts per locus. Those counts are a join and a group-by away in polars; the dict is the only place the function holds the full candidate set in Python.

---

## Verified correct — do not re-flag

**The FUMA algorithm.** All three stages match Nakamura/Watanabe's definitions and the shipped defaults are FUMA's (`lead_pvalue 5e-8`, `candidate_pvalue 0.05`, `clump_r2 0.6`, `lead_r2 0.1`, `window_kb 250`, `merge_distance_bp 250000`):

- **Independent significant SNPs** — greedy clumping in descending LP over variants at `lead_pvalue`, depleting *only* the significant pool at r² ≥ `clump_r2` within the window (`:1138-1182`). Depleting only that pool is correct: FUMA candidates may tag more than one independent SNP, and partners include reference-only variants that were never candidates. The code says so at `:1173-1175`.
- **Candidate membership** — reference-only partners are retained while GWAS-tagged partners above `candidate_pvalue` are dropped (`:1160-1164`). That is FUMA's candidate definition, and `_window_reachable_gwas_annotations`'s docstring (`:1042-1052`) correctly explains why the annotation subset must keep the above-threshold GWAS rows so they are *excluded* rather than misread as reference-only.
- **Lead SNPs** — greedy clumping of the independent heads at r² ≥ `lead_r2` (`:1244-1263`).
- **Loci** — sorted interval merge within `merge_distance_bp` on the same chromosome, with transitive extension through `merged[-1]["end"]` (`:1486-1500`).

**The lead-versus-locus-P invariant holds.** `initial` takes `p.min()` / `lp.max()` across a lead's group (`:1409-1410`) while `l_beta` / `l_se` / `l_rsid` come from the row where `ind_sig_SNP_id == lead_SNP_id`. Those could disagree — but they cannot: `find_lead_snps` consumes heads in descending LP and absorbs only *unassigned* heads, so every member of a lead's group has LP ≤ the lead's. The minimum P and maximum LP are therefore the lead's own. On a distance merge, `_stronger` promotes the whole attribute set together (`:1503-1516`) rather than leaving a mixed row.

**LP is the ranking key wherever P would underflow.** `10^-LP` underflows to 0.0 below ~1e-308 — exactly the range lead SNPs occupy — and the code ranks on LP in `add_canonical_ids` (`:152-159`), both clumping loops (`:1106-1112`, `:1225-1230`) and the merge tie-break (`:1501-1503`), with a comment at each site.

**Allele and identifier handling.** `canonical_variant_id` sorts alleles so A/G and G/A are one variant; `add_canonical_ids` refuses to flip alleles and *raises* on duplicate canonical IDs with conflicting effect-allele orientation (`:218-227`) rather than guessing, keeping the strongest row for genuine duplicates.

**`_window_reachable_gwas_annotations` cannot cross chromosomes** — it is called with `gwas_subset` inside `process_chromosome` (`:1770`), so the `join_asof(strategy="nearest")` compares positions only within one chromosome.

**Contig naming.** `resolve_reference_contig` (`:491-518`) reads the reference's own labels from the tabix index and matches `6`/`chr6` in either direction, turning a genuine naming mismatch into one explicit error. Same guard `sumstat_filter` has; the one `annot_ldblock` lacks.

**Reference-manifest validation is the best input contract in this codebase.** `service.py:229-293` checks format version, genome build, population membership, orientation, both file patterns, both column lists, and — the part that matters scientifically — that the reference was *built* with a window at least as wide as the run requests and r²/MAF floors at least as permissive. A reference computed at r² ≥ 0.2 cannot be silently used for a run that needs r² ≥ 0.1.

**Row-level reference validation.** `_parse_ld_payload` (`:652-673`) checks chromosome, coordinates and r² ∈ [0,1] on every retrieved row and raises rather than filtering silently; `_inventory_lookup` additionally refuses a reference in which an ID maps to multiple canonical variants or vice versa. Index membership is established from the inventory, never inferred from the presence of a pair row (`:437-443`) — which is what preserves FUMA's distinction between an unsupported index and a supported index with no qualifying partners.

**VCF contract.** `service.py:397-419` reads the header, validates the declared build against `modules.ld_clumping.genome_build`, and requests the `<POP>_LDblock` field only when the `region` method is active (`:318`).

**Configuration-source discipline.** `_resolve_standard_options` (`:76-137`) refuses to let a *scientific* value fall back to packaged defaults when the caller already supplied a resolved module configuration — it raises instead. That guard is unusual and worth keeping. It also short-circuits when nothing is missing, so the per-index-SNP call in `get_ld_partners` is cheap.

**Shared primitives.** `common.py` holds `normalise_chromosome`, `chromosome_expression`, `canonical_variant_id` and `exclude_mhc` for both methods, with a docstring stating the reason ("what stops the two methods from disagreeing about what `chr01` means"). `exclude_mhc` returns audited statistics so both callers log identical fields. This is the reuse pattern the rest of the codebase should follow.

**Failure reporting.** Excluded indexes are written to a dedicated TSV with the reason and the observed reference IDs, so an excluded index is explainable rather than merely absent; `PipelineStageError` carries a stable stage/function/context triple throughout; and `_quarantine_changed_outputs` (`service.py:608`) protects partially written outputs on failure.

---

## Suggested order of work

1. **D1** — fix the three example documents to the version-2 six-file-plus-manifest layout. This is the only finding a user hits today, and it currently presents as "the tool ran and found nothing".
2. **I5** — decide whether a reference that covers no chromosome should still exit 0. At minimum, fail when *every* scheduled chromosome was skipped; the per-chromosome `warning` default is reasonable, "warned about all 22" is not.
3. **I2** — vectorise `_inventory_lookup`'s row loop. Measured 11× on 500k rows, and it runs on the largest frame in the module.
4. **I3** — replace `find_lead_snps`'s per-lead frame rescan with the set-based pattern `find_ind_sig_snps` already uses.
5. **I4** — stream tabix output to a temporary file instead of `capture_output=True`.
6. **I1** — revisit the parallelism once 2–4 are done. With the Python loops gone or shrunk, threads may be enough; if not, the frame-passing strategy has to change with it.
7. **D2, D3, O1–O5.**

Nothing here changes a number the module reports. Steps 3–5 change what it costs; steps 1–2 change whether a user can tell that it did not run.
