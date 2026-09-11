# PostGWAS `gcta_cojo` — module review

**Date:** 2026-08-27 · **Scope:** `src/postgwas/modules/gcta_cojo/` (`service.py` 1507, `cli.py` 435, `results.py` 221, `adapters.py` 97, `errors.py` 8), `config/models/modules/gcta_cojo.py`, `config/defaults/modules/gcta_cojo.yaml`, `docs/modules/gcta_cojo/README.md`, and the pipeline wiring through `formatter`.
**Method:** read-only on the tree — no file was modified. Every constructed flag was checked against GCTA-COJO's documented command line and every packaged default against GCTA's own default. The command defect in **I1** was **reproduced** by executing `build_cojo_command` directly. The costs in **I2** are measured on synthetic 5,000,000-variant `.ma` and BIM files, with the validation path split into its three stages. Line numbers refer to the working tree as of 2026-08-27.

**Summary.** The GCTA interface is faithful — correct flags, correct defaults, and the four modes map onto GCTA's four entry points properly, including the detail that `--cojo-p` belongs only to `--cojo-slct`. The input validation is unusually thorough and catches the traps COJO users actually fall into. There are no Critical findings.

Two things are worth acting on: one mode combination emits a **malformed command** (documented as forbidden, never enforced), and the pre-flight validation spends **92% of its time on parsing and index building** for joins that take 8%.

---

## IMPORTANT

### I1 — `joint` mode plus `--cojo-extract` emits two conflicting `--extract` flags

**Where:** `adapters.py:53-56` and `:66-70`. The README (`docs/modules/gcta_cojo/README.md`) states the constraint: "`cojo_extract` … cannot be combined with `joint` mode." Nothing enforces it — not the pydantic model, not `_require_gcta_cojo_arguments`, not `_validate_inputs`, and no test.

Reproduced by calling the builder directly with `mode="joint"`, `joint_snps` and `extract_snps` both set:

```
gcta64 --bfile ref --cojo-file s.ma --maf 0.01 --diff-freq 0.2 --cojo-wind 10000
       --cojo-collinear 0.9 --thread-num 4 --out out
       --extract /tmp/extract.txt --extract /tmp/joint.txt --cojo-joint
                 ^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^
--extract occurrences: 2
```

`--extract` is appended unconditionally whenever `extract_snps` is set, and again by the `joint` branch. Which file GCTA honours is version-dependent and I could not verify it here.

**Severity is contained, and by a good guard rather than by luck.** If GCTA honours the wrong file, the returned SNP set will not match `joint_snps`, and `_validate_modeled_snp_survival` (`service.py:1363-1369`) raises. So this fails loudly rather than producing a wrong joint model — but it fails with a message listing `--cojo-maf`, `--cojo-diff-freq`, `--cojo-chromosome`, extract/exclude lists, alleles and reference genotypes as the likely causes, none of which is the real one. A user would spend a long time on that.

**Fix:** reject the combination in the config model where the README already says it is invalid, and drop the duplicate append in the `joint` branch. A one-line validator and a test.

### I2 — Input validation is 92% parse and index-build, 8% join

**Where:** `_insert_bim` (`service.py:288-336`), `_insert_summary` (`:337-455`), the two SQL checks at `:534-548`.

Both files are read with a per-line Python loop — a regex split per row, five `float()` conversions plus range checks per `.ma` row, and batched inserts into a SQLite table whose `variant_id` is a `PRIMARY KEY`. Measured on 5,000,000 variants in each file:

| stage | time |
|---|---|
| Python parse + per-row checks (both files) | **19.4 s** |
| + SQLite inserts and index maintenance | **~35.6 s** cumulative |
| the two SQL join queries (overlap, allele compatibility) | **4.4 s** |
| **full path as written** | **52.0 s**, peak RSS **50 MB** |

For comparison, the same checks written as vectorised polars over the same files: **10.8 s** — but peak RSS **978 MB**.

**The design choice is defensible and should not be reversed wholesale.** SQLite's 50 MB peak is flat in file size because it streams to disk; the polars version materialises the join, and against a full-WGS BIM rather than a 5M subset that 978 MB becomes many gigabytes. Trading 20× memory for 5× speed is the wrong trade for a module whose inputs are genome-wide by construction.

**What the measurement actually indicts is the row loop, not the storage.** The joins — the part SQLite is genuinely good at — are 8% of the cost. Roughly 31 s of the 52 s is Python parsing and float conversion, and another ~16 s is index maintenance on the inserts. Options that keep the memory profile: read and validate with a streaming vectorised reader and insert only the `(id, allele1, allele2)` triples; or build the two tables without the `PRIMARY KEY` index and add a `UNIQUE` index once after bulk insert (duplicate detection is currently paid per row); or drop the index entirely and detect duplicates with a `GROUP BY` alongside the two joins already being run.

This runs on every invocation, before GCTA starts, and there is no fast path on resume.

### I3 — The strictness of the two reference checks is inconsistent

**Where:** `service.py:534-559`.

- **Allele incompatibility: one variant is fatal.** Any ID present in both files whose alleles do not match in either orientation raises and stops the run.
- **Non-overlap: half the file is a warning.** `minimum_reference_overlap_fraction: 0.50` means up to 50% of `.ma` variants can be absent from the BIM with only a warning appended.

Those are opposite postures towards the same class of problem — a reference that does not match the study. The allele check is the stricter one, and it is also the one most likely to fire spuriously: a single rsID reused for different alleles between releases (common where a multi-allelic site was collapsed) kills the run, where GCTA itself would drop that variant and continue.

There is a second edge worth naming: **a strand-flipped reference makes every variant incompatible.** `.ma` A/G against BIM T/C is the same variant on the other strand, and the check treats it as a mismatch. The resulting error says "incompatible alleles" without mentioning strand, which is the first thing the user should check when the count is very large.

Consider a tolerance on the allele check that mirrors the overlap one — a configurable fraction, fatal above it — and adding "or the reference is on the opposite strand" to the message.

---

## OPTIONAL

**O1 — `--cojo-actual-geno` is not exposed.** GCTA offers it for the case where the GWAS and the LD reference are the same individuals, which changes how COJO computes the variance. It is a niche option and omitting it is reasonable, but the README's reference-requirements section is the place to say so explicitly, since a user with in-sample genotypes would otherwise assume the default is appropriate.

**O2 — `gcta_cojo` is the only module without a wiki page.** `docs/wiki/modules/` has `gcta-gene.md` but no `gcta-cojo.md`; `home.md:102`, `command-reference.md:18` and `README.md:1293` all send users out of the wiki to `docs/modules/gcta_cojo/README.md`. That README is better than several wiki pages, so this is placement, not content — but it makes the module look less finished than it is.

**O3 — `service.py` is 1,507 lines** and holds configuration resolution, four kinds of input validation, SQLite schema management, staging, resume, output-contract checking and two screen renderers. `results.py` and `adapters.py` show the seams are known; the validation block (`_create_validation_database` through `_validate_inputs`, roughly 380 lines) is the natural next extraction.

---

## Verified correct — do not re-flag

**The flag mapping is faithful to GCTA.** `--bfile`, `--cojo-file`, `--maf`, `--diff-freq`, `--cojo-wind`, `--cojo-collinear`, `--thread-num`, `--out`, `--chr`, `--extract`, `--exclude`, `--cojo-gc` with its optional lambda, and the four mode entry points `--cojo-slct`, `--cojo-top-SNPs`, `--cojo-joint`, `--cojo-cond`.

**`--cojo-p` is passed only in `slct` mode** (`adapters.py:62-63`). GCTA ignores it for `cond` and `joint`, and passing it there is a common way to convince yourself a threshold applied when it did not. This module gets it right.

**The packaged defaults are GCTA's own defaults**, not invented ones: `--cojo-p 5e-8`, `--cojo-wind 10000` kb, `--cojo-collinear 0.9`, `--diff-freq 0.2`. The reference sample-size warning at 4,000 individuals matches GCTA's own recommendation for COJO LD references.

**`joint` mode is constructed the way GCTA documents it** — `--extract <joint list>` together with `--cojo-joint` — which is why the duplicate in I1 is a stray append rather than a misunderstanding of the tool.

**Allele compatibility is genuinely validated**, in both orientations, and raises (`service.py:537-548`). The README's claim to validate "`.ma`/BIM ID overlap and allele compatibility" is accurate; I checked specifically because ID-only overlap checks are the usual shortcut here.

**The SNP-list cross-checks are subtle and correct** (`:596-628`):
- conditioning SNPs may not appear in the exclude list;
- **every conditioning SNP must be present in the extract list**, so a narrowing `--extract` cannot silently truncate the conditioning model — this is the non-obvious one, and it is the trap most likely to produce a quietly wrong conditional analysis;
- joint-model SNPs may not appear in the exclude list.

**`_validate_modeled_snp_survival` (`:764-796`) is the module's best guard.** For `cond` and `joint` it requires GCTA's returned SNP set to match the requested list **exactly** — no missing, no unexpected — because a silently dropped conditioning or joint SNP means the model fitted is not the model requested. It is also what contains I1's blast radius.

**Duplicate identifiers are rejected everywhere they matter**: in the `.ma` and the BIM via `PRIMARY KEY` violations with line-numbered errors, and in GCTA's own output via `_read_output_snp_ids` and `normalize_cojo_results`.

**The `.ma` contract is validated per row** — header matching either the formatter's spelling or GCTA's official one, field count, non-empty and distinct alleles, finite numerics, effect-allele frequency strictly inside (0, 1), positive standard error, P within [0, 1], and a minimum sample size.

**Result schemas are validated per mode** (`gcta_cojo.yaml:76-132`) with distinct column sets for `slct`/`top_snps`/`joint` (`bJ`, `bJ_se`, `pJ`) and `cond` (`bC`, `bC_se`, `pC`), and the model validator requires every mode to have a schema and every reported statistic column to be a required column. The `not_estimable_collinearity` status is preserved rather than dropped, so a SNP GCTA could not estimate is visible instead of missing.

**Build and reference population must be declared and are never inferred** (`_validate_build_and_population`, `:170-190`), matching the honest posture taken in the LDSC module for the same unverifiable property.

**Full provenance and resume.** `configuration_digest`, `resolve_completion_resume`, `apply_completion_restart` and `write_completion_manifest` are all used, the registry marks `direct_checkpoint="native"`, outputs are staged and published only after the expected-output contract is checked, and the "no signals" outcome is recorded in the manifest so a legitimate empty result resumes correctly instead of re-running.

**The README's scientific guidance is correct and prominent** — that COJO needs the *unfiltered* summary statistics because it uses them to estimate phenotypic variance, so a region should be restricted with `--cojo-extract` or `--cojo-chromosome` rather than by pre-filtering the `.ma`; and that these modes are not mtCOJO or COJO-SBLUP. Both are exactly the misunderstandings that produce wrong COJO results in practice.

---

## Suggested order of work

1. **I1** — reject `joint` + `extract_snps` in the config model, remove the duplicate append, add the test. Small, and it converts a misleading failure into a clear one.
2. **I3** — decide the intended posture on reference mismatch and make the two checks consistent; add strand to the incompatible-alleles message.
3. **I2** — attack the row loop and the per-row index maintenance, keeping SQLite for the join. Measured: 4.4 s of the 52 s is the part SQLite is doing well.
4. **O2** — give the module a wiki page.
5. **O1, O3.**

Nothing here changes a number the module reports. Step 1 fixes a command it can build but should not.
