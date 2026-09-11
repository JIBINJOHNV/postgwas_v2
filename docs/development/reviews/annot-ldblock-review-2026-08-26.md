# PostGWAS `annot_ldblock` — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/ld_annotation/` (3 files, 150 lines), `config/models/modules/ld_annotation.py`, `config/defaults/modules/ld_annotation.yaml`, `docs/wiki/modules/ld-annotation.md`, and the `annot_ldblock → ld_clump` contract.
**Method:** read-only on the tree — no file was modified by this review. Every behavioural claim about bcftools was verified empirically against **bcftools 1.19 / htslib 1.19** with purpose-built fixtures, not asserted from the manual. Line numbers refer to the working tree as of 2026-08-26.

**Summary.** The science bcftools performs here is correct — coordinate base, block boundaries, boundary-spanning indels and multi-population layering all behave properly (see *Verified correct*). What is wrong is everything around the call: nothing checks that the LD blocks belong to the same genome build or the same chromosome naming as the VCF, and **both failures are completely silent**. Separately, the module performs four full read/write passes and three index builds over a whole-genome VCF where one of each would do.

---

## CRITICAL

### C1 — The genome build of the LD blocks is never checked against the VCF

**Where:** `annot_ldblock.py:54` (and `:20`); `ld_annotation.yaml:2`; `common.py:419-427`.

**What.** `genome_build` is used for exactly one thing — spelling a filename:

```python
bed = Path(ld_dir) / f"{genome_build}_{pop}_ldetect.bed.gz"
```

The CLI default for that value is `GRCh37`, sourced from `modules.ld_annotation.genome_build` in the shipped YAML. Annotate a **GRCh38** harmonised VCF without passing `--genome-build` and you get GRCh37 LDetect blocks applied to GRCh38 coordinates. Every variant receives a label, every label is wrong, and nothing in this module notices — the output is a well-formed VCF full of plausible, incorrect block assignments.

**Why the downstream guards do not save you.**

- `ld_prune_region.py:156-161` raises only when *zero* variants are annotated. Wrong-build annotation produces a full set of labels, so it never fires.
- `ld_clump` does check the declared build (`service.py:401-419`) — but one module later, after the entire VCF has been rewritten, against `modules.ld_clumping.genome_build`, and with an error that blames LD-clumping rather than the annotation that is actually wrong.
- **Standalone** (`postgwas annot_ldblock …`, a documented mode — `docs/wiki/core/running-modules-independently.md:16`) there is no check at any point, ever. The mis-labelled VCF is the deliverable, and anything consuming it downstream (FLAMES, CALDERA, an external tool) inherits wrong loci.

**Fix — six lines, reusing what already exists.** `core/vcf.py:122 validate_vcf_header_contract()` already parses the declared build out of the `##genome_build=<build>` header, and is already imported by `filtering` (`sumstat_filter.py:22`) and `ld_clumping` (`service.py:33`). `annot_ldblock` is the only VCF consumer in this chain that skips it. Copy the ld_clumping pattern verbatim:

```python
header = read_vcf_header(vcf_path, bcftools, ...)
declared_build, contigs = validate_vcf_header_contract(
    header=header,
    genome_build_header=configuration.modules.harmonisation.vcf_processing.genome_build_header,
    supported_genome_builds=list(configuration.resources.genomes),
    required_fields=(),
)
if declared_build != genome_build:
    raise ...   # same message shape as service.py:415-419
```

The harmonisation writer also emits `##contig=<ID=1,length=249250621, assembly=GRCh37>` (`adapters/gwas2vcf/vcf.py:149`), so contig lengths are an independent second check if the `##genome_build` line is ever absent.

### C2 — A chromosome-naming mismatch annotates nothing, silently, with exit status 0

**Where:** `annot_ldblock.py:73-83`.

**Verified.** A BED using `chr1` against a VCF using `1`:

| | observed |
|---|---|
| bcftools exit status | `0` |
| stderr | *empty* |
| output VCF | valid, complete, indexed |
| variants annotated | **0 of 8** |

So `chr`-prefixed blocks against a non-prefixed VCF — or the reverse — produce a well-formed, fully unannotated output that this module reports as success. `docs/wiki/modules/ld-annotation.md:87` already lists "mismatched `chr` naming" as a known problem; it is diagnosed nowhere in the code.

**Where it eventually surfaces.** Only in `ld_clump` when `region` is among the methods (`ld_prune_region.py:156-161`). With `standard`-only clumping the tag is never even queried (`service.py:318`: `if "region" in module.methods and include_region_annotation`), so nothing checks it. Standalone: never.

**Fix.** `validate_vcf_header_contract` already returns the VCF's contig list as its second return value — currently discarded by every caller (`ld_clumping/service.py:401` writes it to `_`). Compare it against the BED's contigs, which `tabix -l <bed.gz>` returns from the index in constant time. An empty intersection raises *before* a multi-GB file is rewritten. Two cheap calls, both already available in the codebase.

---

## IMPORTANT

### I1 — Four full read/write passes and three index builds where one of each is enough

**Where:** `annot_ldblock.py:47`, `:50-89`.

Per run over N populations the module performs: 1 full file copy + N × (decompress + recompress) + N × tabix build. For the default three populations on a multi-GB harmonised VCF that is **four writes and three indexes**.

Two of those three costs are pure waste, both verified:

- **bcftools annotate does not require the *target* VCF to be indexed** — only `-a` benefits from an index. Confirmed on an unindexed target. Only the final index is needed, so **2 of 3 `tabix` builds are dead work** (`:89`).
- **The `cp` at `:47` writes a full copy that the very first `annotate` immediately supersedes.** The first pass can read `vcf_path` directly.

The remaining N passes collapse to one by piping uncompressed BCF between stages. Verified to produce identical annotations:

```bash
bcftools annotate -a EUR.bed.gz -c CHROM,FROM,TO,EUR_LDblock -H '##INFO=…' -Ou in.vcf.gz \
 | bcftools annotate -a AFR.bed.gz -c CHROM,FROM,TO,AFR_LDblock -H '##INFO=…' -Ou \
 | bcftools annotate -a EAS.bed.gz -c CHROM,FROM,TO,EAS_LDblock -H '##INFO=…' \
       -Oz --threads N -o out.vcf.gz
tabix -p vcf out.vcf.gz          # once
```

Net: **1 read + 1 write + 1 index**, roughly a 4× reduction in I/O on the dominant cost of this module.

Serial application is genuinely required — `bcftools annotate` takes one `-a`, and the population BEDs deliberately have *different* intervals, so they cannot be merged into a single `-c` list. Pre-merging the three BEDs into one multi-column annotation resource at resource-prep time would allow a single pass, but that is a new resource to build and version; the pipe gets most of the win with none of it.

One detail for the piped form: `--threads` on `annotate` is *output compression* threads only. Give them to the final `-Oz` stage; the `-Ou` stages get no benefit from them.

### I2 — `os.system(f"cp {current_vcf} {annotated_vcf}")`

**Where:** `annot_ldblock.py:47`.

Five distinct problems in one line:

1. **Exit status ignored.** A copy that fails on a full disk or a permissions error is not detected; the next bcftools failure gets the blame.
2. **Unquoted interpolation into a shell.** Breaks on any path containing a space, and is a shell-injection surface for an operator-supplied path.
3. **Bypasses the module's own `run_cmd`**, imported two lines above, which exists for exactly this and raises with the command and stderr attached.
4. **Does not copy the `.tbi`.** With an empty population list (which the CLI permits — see O5) the output ships without an index.
5. **Not re-entrant.** `run_annot_ldblock_runner` sets `args.vcf = outputs["annotated_vcf"]` (`runners.py:167`). Re-running the step into the same output directory makes source and destination the same path; `cp X X` fails and the failure is discarded.

Under I1 this line disappears entirely, which is the fix.

### I3 — The block-label format is a hard contract with `ld_clump` that nobody declares or validates

**Where:** `ld_prune_region.py:190-206` vs. `common.py:89-99` and `docs/wiki/modules/ld-annotation.md:23`.

The region clumper recovers each block's coordinates by regex **from the label string**:

```python
pl.col(ld_field).str.extract(r"_(\d+)_(\d+)$", 1).alias("START"),
pl.col(ld_field).str.extract(r"_(\d+)_(\d+)$", 2).alias("END"),
```

and raises `"Found %d variants with malformed LD-block coordinates"` if either is null or `END <= START`. So the fourth BED column **must end in `_<start>_<end>`**. The `--ld-region-dir` help text says only "the fourth column is used as the LD-block annotation label"; the module docs say only "annotation-label columns". Neither states the format.

A user-supplied BED with plain identifiers (`chr1_block_17`) therefore annotates perfectly and fails *after the entire VCF has been rewritten*, in a different module, with a message about coordinates the user never knew were encoded in a name.

**Fix.** Test the first few label lines against the same regex during preflight (`zcat bed.gz | head` — cheap), and state the requirement in both the help text and the docs.

### I4 — Dead configuration that promises behaviour the code does not have

- **`include_unassigned`** (`models/modules/ld_annotation.py:9`, `ld_annotation.yaml:4`) is declared, typed and defaulted — and read by nothing in `src/` or `tests/`. Meanwhile `ld_prune_region.py:156` unconditionally drops unassigned variants (`filter(pl.col(ld_field).is_not_null())`) and never counts them. The knob that appears to govern this does nothing, and the behaviour it appears to govern happens silently in another module.
- **`resources.ld_blocks: Path | None`** (`models/resources.py:49`) is likewise declared and never read. It is the obvious default for `--ld-region-dir`, which is currently mandatory on every invocation.

Decide per field: wire it or delete it. Dead config is worse than absent config, because it reads as a documented guarantee.

### I5 — No report of what was annotated

The module returns `{"annotated_vcf": path}` and writes no log. Given C2, **the assignment rate is the one number that distinguishes success from silent total failure.** The docs already concede this (`ld-annotation.md:92-94`).

Note also that two command runners exist in this codebase — `core/processes.run_checked_command` (logger-aware, with `expected_outputs` verification, used by `core/vcf.py` and `ld_clumping`) and the thinner `core/execution/runtime.run_cmd`. This module uses the weaker one and so gets neither logging nor output verification for free.

**Cheapest useful version:** one `bcftools query -f '%INFO/EUR_LDblock\t…\n'` pass over the final VCF, counting non-`.` values per column, logged per population, raising on zero. That is one extra streaming pass — and with I1 applied the module still ends up doing far less total I/O than it does today.

---

## OPTIONAL

**O1 — Dead imports.** `annot_ldblock.py`: `subprocess`, `ThreadPoolExecutor`, `as_completed` are unused (`os` is used only by the line I2 removes). `service.py`: **all six** of `argparse`, `sys`, `RichHelpFormatter`, `validate_path`, `Path`, `require_executable` are unused. The `ThreadPoolExecutor` import is actively misleading — the docstring says "Sequentially annotate", and the passes are order-dependent, so it advertises a concurrency design that cannot exist.

**O2 — `max_memory_gb` is accepted and never used.** `annot_ldblock.py:24`, default `'6G'` — a string, despite the `_gb` name — and `service.py:45` computes `memory_limit(args.memory_gb)` to feed it. Delete both. bcftools annotate's footprint here is bounded by ~1,700 intervals.

**O3 — Non-atomic replace.** `:85-87` unlink then rename; `Path.rename` is already an atomic replace on POSIX, so the explicit unlink only creates a window in which the output does not exist. Between the rename and the `tabix` on the next line the file also carries the *previous* population's stale index. Both disappear under I1.

**O4 — Required-argument validation in three places, naming a flag that does not exist.** `cli.py:71-79` checks manually and exits 1; `service.py:27-33` checks again and raises; neither uses argparse `required=True`. `service.py:28` calls the flag `--ld-dir` — it is `--ld-region-dir`.

**O5 — `--ld-block-populations` accepts duplicates.** `nargs="+"` with `choices` (`common.py:76-88`) permits `EUR EUR`, while `LDAnnotationConfig.unique_populations` forbids it. The duplicate costs a second full annotation pass to write the same tag.

**O6 — The default population list is hardcoded twice**: `annot_ldblock.py:23` and `ld_annotation.yaml:3`.

**O7 — `Path("STUDY_ldblock.vcf.gz").with_suffix(".tmp.vcf.gz")` yields `STUDY_ldblock.vcf.tmp.vcf.gz`.** Harmless and accidental. It is also fixed per dataset, so two concurrent runs sharing a `dataset_id` and output directory collide on the temp file.

**O8 — Provenance vs. byte-reproducibility.** Each pass appends `##bcftools_annotateCommand=…; Date=…`, so outputs are not byte-reproducible across runs. Keep it — provenance is worth more than byte-identity — but `--no-version` is the switch, should anything ever checksum these VCFs.

---

## Verified correct — do not re-flag

All checked against bcftools 1.19 / htslib 1.19 with purpose-built fixtures.

- **BED coordinates are interpreted 0-based half-open**, which is correct for LDetect BEDs. A block `1 100 200` labelled a variant at 1-based 101–200 and left 100 unassigned. **bcftools selects this from the filename suffix**: the identical bytes named `.txt.gz` were read 1-based-inclusive and shifted every block's left edge by 1 bp. The `<build>_<pop>_ldetect.bed.gz` convention is therefore load-bearing for coordinate correctness and deserves a code comment, because it is invisible.
- `-H/--header-line` exists in bcftools 1.19 and accepts the full header string as constructed at `:68`.
- **Re-annotating an already-annotated VCF does not duplicate the INFO definition** — bcftools deduplicates identical header lines.
- **A missing `.tbi` beside the BED does not break annotation** — htslib falls back to reading the file. Likewise, a lexicographically ordered BED (`1, 10, 2, …`) against a karyotypically ordered VCF annotates correctly. Neither is a bug; equally, neither is a guard, so neither is a reason to skip C2's preflight.
- **An indel spanning a block boundary receives exactly one label** — contiguous blocks produce no `Number=1` violation.
- **Serial per-population application is genuinely required** — one `-a` per invocation, and population intervals differ by design.
- The INFO declaration `Number=1,Type=String` matches exactly what the region clumper reads.

---

## Documentation

`docs/wiki/modules/ld-annotation.md` is accurate about what the module does, and unusually honest in its *Limitations* section — it already names the copy-before-validation, the sequential passes and the missing assignment report. Three gaps:

- The required label format `…_<start>_<end>` (I3) appears nowhere, though `ld_clump` hard-depends on it.
- *Input requirements* should state that the file must be named `.bed.gz` for coordinates to be read as BED (see *Verified*), and that build and chromosome naming must match the VCF **because nothing checks them** (C1, C2).
- *QC and logs* advises confirming that "expected variants have block labels" — correct advice with no tooling behind it. If I5 lands, replace it with the logged per-population counts.

---

## Suggested order of work

1. **C1 + C2 preflight.** One `read_vcf_header` + `validate_vcf_header_contract` call (already imported by two sibling modules), compare the declared build to `genome_build`, and intersect the returned contig list with `tabix -l` on each BED. ~15 lines, no new abstraction, and it closes both silent-wrong-answer paths *before* any I/O is spent.
2. **I2 + I1 (partial).** Delete the `cp`, annotate from `vcf_path` on the first pass, move `tabix` out of the loop. Removes one full write and two of three index builds — a net deletion of about five lines.
3. **I3.** Validate the label format in that same preflight; fix the help text and the docs.
4. **I5.** Per-population assignment counts, logged, raising on zero.
5. **I1 (full).** Collapse the N passes into one `-Ou` pipe. The only item needing a new helper — a `run_piped_commands` alongside `run_checked_command` in `core/processes.py`, which harmonisation and formatter would also use. Defer until 1–4 are in.
6. **I4 + O1–O8.** Decide `include_unassigned` and `resources.ld_blocks`, then the dead-code sweep.

Steps 1–4 change what the module can get *wrong*. Step 5 changes what it *costs*.
