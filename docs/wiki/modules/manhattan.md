# Manhattan Plots

## Purpose

The Manhattan command plots GWAS association strength across the genome from a
harmonised GWAS-VCF, with optional consequence and allelic-shift annotations.

## What the analysis does

PostGWAS validates the single VCF sample, constructs the bundled R plotting command,
applies AF and LP display thresholds, optionally marks CSQ or allelic-shift
annotations, and writes a PNG or PDF, the complete R transcript and a point audit.

## When to use it

Use it for visual QC and presentation after harmonisation or filtering. A plot
is exploratory evidence and does not replace formal multiple-testing analysis.

## Input requirements

An indexed single-sample GWAS-VCF with a supported genome-build declaration;
R with compatible optparse, data.table and ggplot2 packages; `bcftools`; and
the required dataset ID/output directory. LP must be numeric Number=1 or
Number=A for biallelic records. Coding highlighting additionally requires
structured CSQ annotations and a working split-vep plugin. An explicit
`--genome-build` or `--pheno` must agree with the VCF declaration.

## Command

```console
postgwas manhattan --vcf PATH [--png PATH | --pdf PATH] [options]
```

## Minimal example

```console
postgwas manhattan \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --png STUDY_manhattan.png \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas manhattan \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --png STUDY_annotated.png \
  --dataset-id STUDY \
  --output-directory results \
  --genome-build GRCh37 \
  --pheno STUDY \
  --min-af 0.01 \
  --min-lp 2 \
  --loglog-pval 10 \
  --csq \
  --width 10 \
  --fontsize 12
```

## Parameters

The canonical `modules.manhattan` YAML supplies 22 autosomes, minimum AF 0,
minimum LP 2, log-log switch at LP 10, cytoband ratio 25, chromosome spacing
20, width 7 inches, and font size 12; height and maximum plot height are unset.
Explicit CLI options override run YAML and packaged YAML. The validated single
VCF sample is used; missing or ambiguous metadata fails instead of being guessed.
`file_format` selects PDF or PNG unless an explicit `--pdf` or `--png` selects
the format/path; these flags are mutually exclusive. `--max-height` is a ceiling
in raw LP units, before display transformation. `axis_label_rows` staggers
chromosome labels; `caption_template` describes filters and display scaling.

## Processing steps

Resolve configuration, validate indexed VCF metadata and the actual R package
stack, query configured fields, validate biallelic point data and assembly
coordinates, render PNG/PDF, and verify the nonempty plot and point audit.
Chromosome offsets are named by assembly chromosome, preserving positions when
intermediate chromosomes are absent. The first headerless query record is retained.

## Outputs

The selected image plus `<dataset>_manhattan.log` and
`<dataset>_manhattan_points.tsv` in the output directory. The default image is
`<dataset>_manhattanplots.<file_format>`. All three relative output patterns are
configurable. The point audit retains raw LP, plotted LP, original chromosome
label/position and cumulative display coordinate.

## QC and logs

Check resolved settings, exact command, sample/build, R/package versions,
retained/excluded counts and point data. Failed or empty queries cannot be
reported as completed plots. Verify labels, caption and association peaks visually.

## Interpretation

LP is `-log10(P)`. The minimum LP option hides weaker associations; it does not
change the underlying VCF. Above threshold t (default 10), plotted LP is
`t * log10(LP) / log10(t)`; t must exceed one. Reference lines and the ceiling use
the same transformation, and the caption identifies log-log spacing. CSQ
highlighting retains unannotated sites as unhighlighted points. Allelic-shift
mode uses a symmetric exact two-sided binomial probability, capped at one.

## Common problems

Missing R packages, wrong phenotype/sample name, VCF without queried fields,
build mismatch, or writing to an unwritable image directory.

## Limitations

Only chromosomes 1–22 and X are supported by the bundled assembly lengths.
The configured `strip_prefix` label policy removes a leading `chr` for plotting,
records affected rows, and preserves the original label in the audit; `exact`
disallows this normalization. Unsupported M/MT and multiallelic plotted records
fail. Allelic-shift mode cannot be combined with the AF filter. This module
produces a Manhattan plot, not a QQ plot.

## Scientific references

- [GWAS-VCF specification](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
- [bcftools query](https://samtools.github.io/bcftools/bcftools.html#query)
- [split-vep keep-sites behavior](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html)
- [UCSC GRCh37 chromosome lengths](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes)
