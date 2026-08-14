# Manhattan Plots

## Purpose

The Manhattan command plots GWAS association strength across the genome from a
harmonised GWAS-VCF, with optional consequence and allelic-shift annotations.

## What the analysis does

PostGWAS selects a phenotype/sample, constructs the bundled R plotting command,
applies AF and LP display thresholds, optionally marks CSQ or allelic-shift
annotations, and writes a PNG or PDF plus the complete R transcript.

## When to use it

Use it for visual QC and presentation after harmonisation or filtering. A plot
is exploratory evidence and does not replace formal multiple-testing analysis.

## Input requirements

A GWAS-VCF; matching genome-build declaration; R with the packages required by
the bundled script; `bcftools` for phenotype auto-detection; and an output path
or dataset/output directory for the default PDF.

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
CLI options override those values. When no phenotype is given, the first VCF
sample is used; if detection fails, PostGWAS passes `unknown_sample`. When
neither image option is given, it writes
`<output>/<dataset>_manhattanplots.pdf`.

## Processing steps

Resolve phenotype, choose PNG or PDF, assemble arguments for
`resources/assoc_plot.R`, create the output directory, write the command header,
and run `Rscript` with stdout/stderr captured.

## Outputs

The selected image and a timestamped
`<image-stem>.<phenotype>.assocplot_<timestamp>.log` beside it.

## QC and logs

Check the logged command, phenotype, build, thresholds, R output, and image
existence. Verify chromosome labels and expected association peaks visually.

## Interpretation

LP is `-log10(P)`. The minimum LP option hides weaker associations; it does not
change the underlying VCF. CSQ highlighting depends on compatible annotation.

## Common problems

Missing R packages, wrong phenotype/sample name, VCF without queried fields,
build mismatch, or writing to an unwritable image directory.

## Limitations

If both `--png` and `--pdf` are supplied, PNG wins; choose only one image
format.

## Scientific references

- [GWAS-VCF specification](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
