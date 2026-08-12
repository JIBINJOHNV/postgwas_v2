# QC Summary

## Purpose

The standalone QC command calculates summary metrics for an existing GWAS-VCF
without producing a replacement filtered VCF.

## What the analysis does

It runs bcftools statistics, extracts essential record/type metrics, and counts
missing study/reference AF values, comparable AF pairs, and pairs whose absolute
difference exceeds the selected cutoff.

## When to use it

Use it for a quick summary of a harmonised or filtered VCF. Use harmonisation's
integrated assessment when you need the full raw-versus-virtual-filter report.

## Input requirements

A GWAS-VCF, dataset ID, output directory, `bcftools`, and the correct
`INFO/<reference AF>` tag name.

## Command

The pipeline target is `qc_summary` and the standalone command is `postgwas qc`.

```console
postgwas qc --vcf PATH [options]
```

The standalone command does not currently expose a `--bcftools` option even
though the QC service requires one, so it stops before reading the VCF. Until
that is corrected, run the QC summary through the pipeline and pass `--bcftools`
explicitly.

## Minimal example

```console
postgwas pipeline \
  --modules qc_summary \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --bcftools bcftools
```

## Full example

```console
postgwas pipeline \
  --modules qc_summary \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --bcftools bcftools \
  --reference-af-column EUR \
  --maximum-af-difference 0.2 \
  --threads 8
```

## Parameters

The command defaults to reference tag `EUR` and maximum absolute AF difference
0.2. Compute options are shared with other commands. This module does not read
`modules.qc_summary` from a run configuration; the values above come from the
command line.

## Processing steps

PostGWAS runs `bcftools stats`, parses essential metrics, evaluates AF
missingness/concordance, transposes the metrics into a two-column report, and
writes it as TSV.

## Outputs

`<output>/<dataset>_qc_summary.tsv` contains metric names and `raw_variants`
values. The bcftools statistics file is written as `<vcf>.stats`, that is, beside
the input VCF rather than in the output directory. The input directory therefore
has to be writable.

## QC and logs

Review total/type statistics, `external_af_missing`, `study_af_missing`,
`af_comparable`, and `af_difference_above_cutoff`. A metric can be unavailable
when the VCF lacks a required field.

## Interpretation

AF-difference counts are meaningful only when the study AF and selected
reference AF are allele-aligned and population/build compatible.

## Common problems

Wrong INFO tag, missing FORMAT/AF, incompatible AF population, unindexed or
invalid VCF, or missing bcftools.

## Limitations

This command is a compact summary. It does not materialize a filtered VCF,
apply the filtering module, or reproduce the full harmonisation QC assessment.
The standalone `postgwas qc` entry point is currently unusable because it cannot
accept the `bcftools` path the service requires.

## Scientific references

- [bcftools statistics documentation](https://samtools.github.io/bcftools/bcftools.html#stats)
- [GWAS-VCF specification](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
