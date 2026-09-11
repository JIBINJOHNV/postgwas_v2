# LD Annotation

## Purpose

LD annotation adds population-specific, approximately independent LD-block
labels to a harmonised GWAS-VCF.

## What the analysis does

For each requested population, PostGWAS annotates intervals from a compressed
BED file into a new VCF INFO field named `<POP>_LDblock` by default.
Population-specific bcftools processes remain sequentially dependent, but they
exchange uncompressed BCF through pipes. Only the final result is compressed
and indexed.

## When to use it

Use it before region-based LD pruning and when downstream logic needs the
Berisa–Pickrell LD-block label carried with each variant.

## Input requirements

A BGZF-compressed, indexed, single-sample harmonised GWAS-VCF with exactly one
supported `##genome_build` declaration, `bcftools`, and one BED file per
population named exactly `<BUILD>_<POPULATION>_ldetect.bed.gz`. A lone VCF
sample ID that differs from the run dataset ID is used with a warning;
zero-sample and multi-sample VCFs fail before annotation. PostGWAS reads `BUILD`
from the VCF header and `POPULATION` from `--ld-block-populations`. For example,
a GRCh38 VCF with EUR and AFR selected requires
`GRCh38_EUR_ldetect.bed.gz` and `GRCh38_AFR_ldetect.bed.gz`. The exact
`.bed.gz` suffix preserves BED coordinate interpretation by bcftools; a tabix
index is supported but is not required for the BED files. The input VCF
requires a valid TBI or CSI so its record count can be read without another
full VCF scan and compared with the annotated output count.

The first four tab-separated BED columns are required:

| Column | Meaning |
|---|---|
| `CHROM` | Exact input-VCF contig name, such as `1` or `chr1`. |
| `START` | Zero-based block start; included. |
| `END` | Zero-based block end; excluded. |
| `BLOCK_LABEL` | Nonempty value written to `<POP>_LDblock`. For region clumping, it must end with `_<START>_<END>`, using the same values as the BED coordinates. |

Additional trailing columns are allowed but ignored. For example, in
`1<TAB>16103<TAB>2047216<TAB>EAS-1_16103_2047216<TAB>LDBLOCK_1`, the first
four columns are used and `LDBLOCK_1` is ignored.

> **Warning:** Every BED file must have been generated for the same genome build
> declared inside the input GWAS-VCF. BED coordinates do not contain enough
> information for PostGWAS to infer or independently verify their genome build.
> Chromosome labels must also match the VCF exactly: `1` matches `1`, and
> `chr1` matches `chr1`, but `1` does not match `chr1`.

## Command

```console
postgwas annot_ldblock --vcf PATH --ld-region-dir PATH [options]
```

`--vcf`, `--ld-region-dir`, `--dataset-id`, and `--output-directory` can
instead be supplied by the canonical LD-annotation run configuration. An
explicit command-line value overrides the corresponding YAML value.

## Minimal example

```console
postgwas annot_ldblock \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas annot_ldblock \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR AFR EAS \
  --dataset-id STUDY \
  --output-directory results \
  --threads 4
```

## Parameters

Canonical `modules.ld_annotation.inputs.ld_region_dir` supplies the shared
directory containing build- and population-specific BED files. For example:

```yaml
modules:
  ld_annotation:
    inputs:
      vcf: data/STUDY.vcf.gz
      ld_region_dir: reference/ld_blocks
      dataset_id: STUDY
    output_directory: results
    populations: [EUR, AFR, EAS]
```

This configuration can be run with:

```console
postgwas annot_ldblock --run-config ld_annotation.yaml
```

The packaged canonical YAML supplies populations EUR, AFR, and EAS. Its
`bed_filename_template`, `info_field_template`,
`info_description_template`, `output_filename_template`, and
`summary_filename_template` settings supply the validated resource, VCF-tag,
provenance-description, annotated-VCF, and QC-summary output patterns. The
`html_report_filename_template` and `canonical_log_filename_template` settings
control the detailed HTML-report and audit-log names.
When annotation and region clumping are selected in one pipeline, the
annotation INFO template must match `modules.ld_clumping.vcf_fields.ld_block`;
configuration validation rejects a mismatch before analysis starts.
The genome build is not a module setting or CLI option: it is read from the
validated input VCF header. `--threads` is applied to final BGZF compression;
the population annotation stages are linked in their required order.

The former `resources.populations.<POP>.ld_blocks` field was unused and has
been removed. Replace any custom occurrence with
`modules.ld_annotation.inputs.ld_region_dir` pointing to the shared directory,
not to a single population BED file; the validated filename template selects
the correct build and population inside that directory.

## Processing steps

Validate the selected population list, read the VCF build, contig declarations,
and indexed record count, locate and validate every selected population BED
before creating output, start one bcftools annotation stage per
population, pipe intermediate stages as uncompressed BCF (`-Ou`), and write the
final stage once as compressed VCF (`-Oz`) with a TBI index
(`--write-index=tbi`). PostGWAS checks every process, validates the temporary
VCF metadata and index, then performs one constant-memory `bcftools query` scan
to calculate annotation coverage and LD-block use. PostGWAS renders the CSV and
self-contained HTML report from the same in-memory coverage evidence without
another VCF scan. The VCF, TBI, CSV, and HTML report are published as one
validated artifact set. A failure leaves any existing output set unchanged and
removes incomplete
temporary files. When detailed progress is enabled, the terminal and durable
screen log show these six validated stages without estimating unobservable
per-variant progress inside bcftools.
Every completed stage prints a concise outcome using values already obtained by
that stage: input build, contigs, and total variants; requested populations and
per-population reference inventories; created INFO fields; output-VCF and index
validation, including the invariant that annotation preserves the input record
count; report preparation; and artifact-publication count. The final terminal
summary is the single detailed presentation of overall and per-population
coverage and all saved output names. Full paths, stage durations, bcftools
commands, detailed coverage evidence, and failure details are written to the
canonical log.

## Outputs

By default:

- `<output>/<dataset>_ldblock.vcf.gz` and its `.tbi` index;
- `<output>/<dataset>_ldblock_summary.csv`;
- `<output>/<dataset>_ldblock_report.html`;
- `<output>/<dataset>_ldblock.log`.

INFO tags default to `<POP>_LDblock` and contain the fourth BED column. The CSV
contains an overall row, one row per selected population, and one row for every
contig with variants unassigned in all selected populations. Population rows
report total, annotated, and unassigned variants, annotation percentage, BED
block count, blocks represented by at least one variant, and empty blocks.
The self-contained HTML report presents the same validated overall and
population-level counts, the contig distribution of variants unassigned in
every population, BED-reference paths, run metadata, output links, and
interpretation guidance. It separates unassigned variants on contigs present
in at least one requested BED from those on contigs absent from every requested
BED and reports how many requested BED files contain each listed contig. BED
availability is a contig-level inventory and does not imply that a particular
variant overlaps an LD-block interval. The report also states explicitly that
annotation is
non-filtering: unmatched variants remain in the output VCF with empty requested
LD-block INFO fields. It does not recalculate coverage independently.

## QC and logs

PostGWAS confirms every requested INFO header has the expected String/Number=1
definition, the genome build and contig declarations are preserved, and the
TBI and streaming QC scan return the same record count. The terminal prints a
MAGMA-style summary of total variants, variants unassigned in every selected
population, unassigned variants on covered BED contigs, uncovered contigs, and
per-population annotation and block coverage. For multiple populations,
"assigned in any population" is the union, "assigned in every population" is
the intersection, and population counts must not be added because the same
variant can carry several population INFO tags. The same findings are preserved
in the CSV, HTML report, and canonical log. The HTML report describes the
supplied resources but does not infer whether their population labels or
coordinate provenance are scientifically correct; that remains the user's
responsibility.

## Interpretation

LD blocks are population- and build-specific approximations. A block label is
not a causal locus and should not be transferred between populations.

## Common problems

Missing or unsupported VCF build metadata, a missing or invalid input VCF
index, wrong filename pattern, mismatched `chr` naming, BED/build mismatch,
missing population file, or unavailable bcftools.

## Limitations

The current engine processes populations sequentially. The sequential
population stages are required because each stage consumes the VCF emitted by
the prior stage; their intermediate files have been eliminated. Summary
calculation requires one additional read-only streaming pass over the final
temporary VCF but does not create a genome-scale intermediate table.

## Scientific references

- [Berisa and Pickrell 2016, approximately independent LD blocks](https://doi.org/10.1093/bioinformatics/btv546)
- [bcftools annotate: BED coordinates and exact sequence-name matching](https://samtools.github.io/bcftools/bcftools#annotate)
- [bcftools performance guidance: use uncompressed BCF between piped subcommands](https://samtools.github.io/bcftools/howtos/scaling.html)
- [bcftools common options: final-output `--write-index`](https://samtools.github.io/bcftools/bcftools#common_options)
