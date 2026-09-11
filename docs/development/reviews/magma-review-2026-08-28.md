# PostGWAS `magma` — module review

**Date:** 2026-08-28 · **Scope:** `src/postgwas/modules/magma/` (`analysis.py` 2137, `service.py` 551, `annotations.py` 569, `cli.py` 319, `reference.py` 100), `config/models/modules/magma.py`, `config/defaults/modules/magma.yaml`, `docs/wiki/modules/magma.md`, `docs/modules/magma.md`, `tests/test_magma.py`. **`magmacovar` and the single-cell MAGMA methods are out of scope**, as requested.
**Method:** read-only on the tree — no file was modified. Every constructed flag was checked against MAGMA's documented command line, every packaged default against MAGMA's own default, and the four multiple-testing corrections against their statistical definitions. The pandas separator behaviour in **I3** was executed, not assumed. Line numbers refer to the working tree as of 2026-08-28.

**Summary.** The MAGMA interface is correct in the places that are easy to get wrong and invisible when wrong: the annotation window operands are in MAGMA's documented order, batching uses MAGMA's own `--batch`/`--merge` rather than a home-rolled gene split, `duplicate=error` refuses silent duplicate handling, and the analysis is seeded. The multiple-testing code is textbook-correct, including a numerically stable Šidák. There are no Critical findings.

**Implementation status (2026-08-31).** The findings below describe the
2026-08-28 snapshot and are retained as review history. The current tree now
has an explicit, build-specific MHC policy and audited SNP/gene exclusions. It
also gives gene-set results one YAML-configured primary global correction
(`fdr_bh` by default), writes the primary method/adjusted p-value/decision into
every result row, and uses only that correction for terminal significance
declarations. Other corrections remain explicitly secondary. Cross-mapping
multiplicity remains a documented study-design responsibility when several
mappings are selected. The delimiter mismatch is resolved: literal
one-character separators and `\s+` retain pandas' C engine, while other valid
multi-character regexes use the Python engine. The pandas-to-Polars result
conversion discussed below remains a separate optional cleanup.

The three findings are about what the module does *not* say: it has **no MHC policy** while every sibling module has one; it emits **up to eight corrected p-values per gene set with no stated primary**; and it mixes pandas and polars in a way that turns one configuration field into a trap.

---

## IMPORTANT

### I1 — No MHC policy anywhere, and MHC inclusion silently depends on another module's setting

**Where:** nothing matches `mhc` in `src/postgwas/modules/magma/`, `config/defaults/modules/magma.yaml`, or `docs/wiki/modules/magma.md`.

Every sibling module that touches genomic coordinates makes this explicit and build-specific:

| module | MHC handling |
|---|---|
| `qc_summary` | `remove_mhc: true`, GRCh37/GRCh38 intervals, counted as a named rule |
| `sumstat_filter` | `remove_mhc: true`, BED exclusion with contig-name matching |
| `ld_clumping` | `remove_mhc: true`, shared `common.exclude_mhc` with audited statistics |
| **`magma`** | **absent** |

This matters more than an omitted default. The extended MHC carries long-range LD that violates the assumptions behind MAGMA's gene-level statistic, and MHC genes routinely dominate the top of a MAGMA gene list. Whether they are present here is decided **elsewhere**: MAGMA consumes the formatter's `snp_loc`/`p_values` files, the formatter has no MHC keys at all, and `sumstat_filter` — the step that would have removed them — ships `enabled: false`. `ld_clumping`'s MHC exclusion is an in-memory filter on its own frame and never reaches the formatter's output.

So with the shipped defaults MAGMA sees the MHC, nothing in its configuration or its results says so, and a user who enables `sumstat_filter` changes the MAGMA gene list without touching anything MAGMA-related.

The right fix is not necessarily to exclude it — including the MHC is what running MAGMA out of the box does, and it is a defensible default. The fix is to make the decision visible and local: an explicit `remove_mhc` with the same build-specific intervals the three sibling modules already share, defaulted to whichever behaviour the project intends, and recorded in the provenance columns so a reader of `03_results/` can tell.

### I2 — Up to eight corrected p-values per gene set, with no stated primary, and nothing correcting across mappings

**Where:** `analysis.py:768-791`; `magma.yaml:98-119`; `docs/wiki/modules/magma.md:136-140`.

The correction families overlap by construction. `go_kegg_reactome` matches `^(?:GOBP_|GOCC_|GOMF_|KEGG_|REACTOME_)`, and `go` matches a subset of the same. A single GOBP gene set therefore receives:

- four global corrections (`bonferroni`, `sidak`, `holm`, `fdr_bh`), plus
- `go_kegg_reactome` Bonferroni and BH, plus
- `go` Bonferroni and BH

— eight adjusted values for one test, all written side by side into `03_results/<mapping>/…_gene_sets_corrected.tsv`. The Interpretation section says only that competitive results "require multiple-testing interpretation"; it never names which column is the intended primary.

Offering several corrections is reasonable for a pipeline whose user is an analyst. Offering eight without saying which one the analysis is reported on is an invitation to pick the significant one after the fact, which is exactly the error multiple-testing correction exists to prevent.

**The second half is larger.** The resource preparer ships **88 mapping definitions** (`docs/wiki/modules/magma.md:104`) and `mapping.selected` is a list, so running positional alongside eMAGMA, H-MAGMA and nMAGMA is an expected use. Corrections are computed **per mapping** — correctly, since the output layout keys on `{analysis}` — but nothing corrects, or even counts, the multiplicity *across* mappings. A user running four mappings has four independent chances at significance for the same gene set and no column reflects that.

Name a primary correction in the docs, and say plainly whether the mappings are alternatives to compare (in which case `04_comparisons/` is the intended reading) or a family that needs its own accounting.

### I3 — pandas, polars and numpy in one module, and a configuration field that only accepts one value

**Where:** `analysis.py:17-19` (all three imported), `_read_table` at `:216-234`, return annotations at `:737` and `:818`.

```python
def _read_table(...) -> pd.DataFrame:
    return pd.read_csv(file_path, sep=delimiter_pattern, comment=comment_prefix, engine="c", ...)
```

Executed against pandas to confirm the behaviour:

```
sep='\s+'   with engine='c'  -> OK          (pandas special-cases this one pattern)
sep='[|;]'  with engine='c'  -> ValueError: the 'c' engine does not support regex separators
```

Four configuration fields are presented to users as regular expressions — `input.table_delimiter_pattern`, `chrom_magma_mapping.mapping_delimiter_pattern`, `.location_delimiter_pattern`, `gene_sets.membership_delimiter_pattern` — and every one of them ships as `\s+`, the single regex the C engine accepts. Configure any other pattern and the run dies with a raw pandas `ValueError` that names neither the field nor the file. It is latent today, but the field's name promises something the reader cannot have.

Separately, `correct_gene_set_p_values` and `correct_gene_p_values` are both annotated `-> pl.DataFrame` while operating on a pandas frame throughout and converting only at `return pl.DataFrame(frame.to_dict(orient="list"))`. The annotations are accurate about the return and misleading about the body, and the round-trip through a Python dict materialises the whole table twice. For a gene-set result this is small; for the gene table it is every gene in the annotation.

---

## OPTIONAL

**O1 — asymmetric collision check.** `correct_gene_p_values` raises if a configured correction column already exists in MAGMA's output (`:849-854`); `correct_gene_set_p_values` performs no such check and would silently overwrite (`:770`). Same hazard, one guard.

**O2 — an underflowed MAGMA p-value of 0 becomes a corrected 0.** `_validate_p_values(..., allow_zero=True)` admits zeros, and every correction maps 0 to 0. That is mathematically right and is MAGMA's limitation rather than this module's, but the output documentation should say that a corrected `0` means "at or below MAGMA's reporting floor", not "exactly zero".

**O3 — `analysis.py` is 2,137 lines** holding command construction, batch planning, both correction routines, gene-set preparation, annotation validation and report writing. `annotations.py` and `reference.py` show the seams are understood; the two correction functions and their table I/O are the cleanest next extraction, and they are the part most likely to be reused by `magmacovar`.

---

## Verified correct — do not re-flag

**The annotation window operands are in MAGMA's order.** `--annotate window=<upstream>,<downstream> --snp-loc … --gene-loc … --out …` (`analysis.py:1546-1553`), with `upstream` resolved from `annotation_window_upstream_kb` falling back to `gene_window_upstream_kb`, and downstream likewise (`:1535-1544`). Reversing these two operands would systematically mis-assign SNPs to genes and would be invisible in every output file; it is right, including the per-mapping override.

**The gene-analysis command is complete and conservative** (`:1209-1221`):

```
--bfile <ref> --gene-annot <annot> --pval <file> use=SNP,P ncol=N_COL duplicate=error
--gene-model snp-wise=mean --seed <execution.random_seed>
```

`duplicate=error` refuses MAGMA's silent duplicate handling; `--seed` from the shared execution seed makes the snp-wise model reproducible run to run. Both are choices, and both are the careful ones.

**The packaged defaults are deliberate and identifiable.** `gene_model: snp-wise=mean` is MAGMA's own default. The `35 kb` upstream / `10 kb` downstream window is *not* MAGMA's default (which is 0/0) — it is the FUMA convention, and it is exposed on the CLI as `--window-upstream` / `--window-downstream` rather than buried.

**Batching uses MAGMA's native mechanism.** `--batch <n> <total>` per worker followed by `--merge <prefix>` (`:1240-1275`). This is the detail most likely to be got wrong by hand: splitting genes across independent runs and concatenating the outputs would produce a `.genes.raw` whose gene-gene correlation structure is not what the gene-set step requires. Using `--batch`/`--merge` avoids that entirely. Parallelism is across MAGMA subprocesses, so — unlike `ld_clumping` — the GIL is not a factor, and a failed batch cancels the rest.

**`_batch_plan` is sensibly bounded** (`:716-729`): workers are the minimum of the thread budget and `memory_gb // memory_per_process_gb`, and batches are additionally capped by `genes // minimum_genes_per_batch`, so a small annotation does not spawn many near-empty MAGMA processes.

**`adjust_p_values` (`core/statistics.py:58-86`) is textbook-correct for all four methods**, and shared rather than duplicated:
- Bonferroni `min(np, 1)`;
- **Šidák computed as `-expm1(n · log1p(-p))`** — the numerically stable form; the naive `1 − (1−p)^n` loses precision at exactly the small p-values that matter;
- Holm as a step-down with `maximum.accumulate((n − i)·p₍ᵢ₎)`;
- Benjamini–Hochberg as a step-up with a reversed `minimum.accumulate(p₍ᵢ₎·n/i)`;
- a stable `mergesort` so ties are deterministic, and results scattered back to input order;
- inputs validated as one-dimensional, finite and within [0, 1] before any of it.

**The correction is auditable and the test count is trustworthy.** Duplicate gene IDs and duplicate gene-set full names are rejected *before* correcting (`:762-765`, `:838-841`), so the multiplier is right; each family's test count is logged together with the pattern that selected it (`:772-791`); and family corrections write `NaN` outside the family rather than a misleading 0 or 1.

**Corrections are per mapping, not pooled.** The output layout keys `corrected_genes` and `corrected_gene_sets` on `{analysis}`, and the call sites sit inside the per-mapping loop (`:1888`, `:2036`).

**The gene-set command is correct**: `--gene-results <genes.raw> --set-annot <prepared sets> --out <prefix>` (`:1982-1984`).

**The output layout is the clearest in the codebase** — `00_run_metadata/`, `01_inputs/`, `02_intermediates/<mapping>/`, `03_results/<mapping>/`, `04_comparisons/`, `05_logs/` — and every mapping's provenance carries `mapping_method`, `biological_context`, `gene_id_type`, `annotation_source`, `annotation_version` and `result_statistic_type`. That last field is what lets `pipeline/runners.py:_validated_magma_gene_result` refuse to hand an uncalibrated mapping to a downstream consumer, which is a real safeguard built on this module's honesty about what each mapping produces.

---

## Documentation

`docs/wiki/modules/magma.md` is accurate on inputs, resume behaviour and failure modes. Three additions:

- **State the MHC position** (I1) — whether the MHC is included, and that with the shipped defaults it is, because no module upstream of MAGMA removes it unless `sumstat_filter` is enabled.
- **Name the primary correction** (I2), and say whether multiple selected mappings are alternatives to compare or a family requiring further accounting.
- **Explain a corrected p-value of 0** (O2) as MAGMA's reporting floor rather than an exact zero.

---

## Suggested order of work

1. **I1** — add an explicit MHC policy with the build-specific intervals the three sibling modules already share, and record it in provenance. Whatever the chosen default, the point is that the decision stops being made by a different module's config.
2. **I2** — name the primary correction in the docs; decide and document the cross-mapping position.
3. **I3** — either validate the delimiter fields against what `engine="c"` accepts, or read with a reader that honours the advertised regex; drop one of pandas/polars from the correction path.
4. **O1–O3.**

Nothing here changes a number the module computes. I1 changes which numbers a user gets without knowing why.
