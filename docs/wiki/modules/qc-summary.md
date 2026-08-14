# QC Summary

## Purpose

The QC module assesses an existing genotype-free GWAS-VCF and writes raw,
rule-level, and combined virtual-QC metrics. It does not alter the input or
create a filtered VCF.

## What the analysis does

PostGWAS evaluates every configured QC rule independently against every raw VCF
record, then reports the conjunction of those rules as a virtual passing subset.
It measures P-value, allele-frequency, imputation-quality, variant-type,
palindromic-SNP, MHC, sample-size, and study/reference frequency evidence without
silently discarding or rewriting variants.

Harmonisation calls the same QC-owned assessment internally, so direct and
integrated reports use one scientific policy and one implementation. Reported
rule counts may overlap; the final accounting separately records unique failures
and records matching multiple rules.

PostGWAS deliberately does not run `bcftools stats` here. That command normally
derives allele-frequency bins and singleton counts from `INFO/AC` and `INFO/AN`,
or from `FORMAT/GT`. Genotype-free GWAS-VCFs normally carry summary-statistic
`FORMAT/AF` instead, so those default bins would not describe the study
frequency. This module reports only metrics derived from explicitly configured
GWAS-VCF fields.

## When to use it

Use QC summary after harmonisation to assess data quality, understand the effect
of candidate filtering rules, or compare raw and virtually passing variants.
Use `postgwas sumstat_filter` when the analysis requires a materialized filtered
VCF. Do not use QC summary to repair an incorrect genome build, allele
orientation, reference population, or sample mapping.

## Input requirements

- A non-empty genotype-free GWAS-VCF with the configured FORMAT and INFO fields.
- A compatible genome build and build-specific MHC policy.
- A dataset ID and output directory.
- `bcftools`, including a version response and access to the input VCF.
- A study sample matching the dataset ID, or an input containing exactly one
  sample.
- A reference-frequency INFO tag whose allele orientation and population are
  compatible with the study frequency when concordance is assessed.

## Command

The standalone command is `postgwas qc`; the pipeline target is `qc_summary`.
Inspect the current direct options with:

```console
postgwas qc --help
```

## Minimal example

```console
postgwas qc \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

Export and edit the canonical YAML, then supply only the run-specific inputs on
the command line:

```console
postgwas config export \
  --module qc_summary \
  --style full \
  --output qc_summary.yaml

postgwas qc \
  --run-config qc_summary.yaml \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --reference-af-column EUR \
  --dataset-id STUDY \
  --output-directory results
```

Explicit CLI values override the corresponding resolved YAML value. Omitted
CLI options do not replace YAML settings.

## Parameters

`modules.qc_summary` is the canonical schema-validated configuration. Export it
instead of copying defaults from an older run. Its main groups are:

- `inputs`, `target_build`, `reference_af_column`, and `output_directory`;
- `rules`, including P-value, AF, INFO, frequency-difference, indel,
  palindromic, MHC, and effective-sample-size policies;
- `vcf_fields`, which defines the exact bcftools queries;
- `table`, which defines temporary-table parsing; and
- `output_layout`, which defines every report, log, configuration, temporary,
  and completion-manifest path.

Missing-value actions are explicit for P value, AF, and INFO. MHC intervals are
one-based and inclusive in the configuration. The packaged Genome Reference
Consortium regions are GRCh37 `chr6:28,477,797-33,448,354` and GRCh38
`chr6:28,510,120-33,480,577`; change the appropriate build entry only when the
analysis protocol justifies a different interval.

Former harmonisation-only QC/filter keys are not a second source of truth.
Equivalent behavior now belongs under `modules.qc_summary`, including
`rules.maf_min`, `rules.info_min`, `rules.info_max`,
`rules.maximum_af_difference`, the missing-value actions, palindromic policy,
MHC regions, and `rules.sample_size_outlier_standard_deviations`.

## Processing steps

1. Resolve and schema-validate the global and QC-summary configuration.
2. Validate the VCF, dataset/sample selection, output paths, configured fields,
   reference-frequency tag, and bcftools executable/version.
3. Validate or restart from the global checksum-backed completion boundary.
4. Run one full-record `bcftools query` into an isolated temporary table.
5. Use streaming Polars aggregation to calculate raw metrics, every independent
   rule, combined pass/fail status, and virtual-subset metrics in one pass.
6. Reconcile raw, passing, uniquely failing, and multiply failing counts.
7. Publish the reports, resolved configuration, completion provenance, log
   records, and concise terminal summary.

The temporary table is removed whether assessment succeeds or fails. For raw
and virtual-QC records, metrics include record and SNP/non-SNP counts,
transition/transversion ratio, AF and INFO missingness, comparable and
discordant frequency pairs, and effective-sample-size availability, range,
mean, sample standard deviation, outlier threshold, and upper-tail count.

## Outputs

With the packaged output layout, the requested output directory contains:

| Output pattern | Description | Important content | Interpretation |
|---|---|---|---|
| `qc_summary/<dataset>_<build>_vcf_qc_metrics.tsv` | Raw and virtual-subset metrics | Counts, missingness, Ti/Tv, AF concordance, effective sample size | Descriptive QC summaries, not filtered data |
| `qc_summary/<dataset>_<build>_vcf_qc_rule_results.tsv` | Per-rule audit | Independent rule counts and combined accounting | Rules can overlap |
| `qc_summary/<dataset>_<build>_qc_assessment.json` | Machine-readable assessment | Complete metrics and balanced accounting | Canonical reusable result contract |
| `qc_summary/<dataset>_<build>_qc_resolved.yaml` | Resolved run configuration | Effective global, resource, and module values | Reproduction metadata |
| `qc_summary/<dataset>_<build>_qc.log` | Canonical service log | Inputs, parameters, decisions, outputs, and final status | Primary execution audit |
| `qc_summary/<dataset>_<build>_qc_complete.yaml` | Completion manifest | Configuration digest and input/output fingerprints | Evidence required for validated resume |

No `.stats` file is written beside the input VCF.

## QC and logs

The canonical log records the resolved rules and fields, VCF metadata,
bcftools path and version, command execution, row accounting, report paths,
warnings, restart decisions, and final completion status. The completion
manifest fingerprints the input VCF and the three scientific assessment reports.

Global resume is enabled by default. A checksum-valid completed result is reused.
A real partial result or changed parameter/input warns and restarts from the
safe QC-summary boundary, replacing only checksum-matching PostGWAS-owned files.
An externally modified output is preserved and refused instead of being deleted
automatically. `--overwrite` remains the explicit forced replacement.

## Interpretation

Reference-frequency concordance is meaningful only when study and reference
frequencies use the same effect/alternate-allele orientation, genome build, and
compatible population. “QC-passed” describes the combined virtual assessment,
not a new VCF and not proof that the study is scientifically suitable for every
downstream method.

## Common problems

| Problem or error | Likely cause | How to check | Solution |
|---|---|---|---|
| Required VCF field is absent | The input schema or configured `vcf_fields` does not match | Read preflight and bcftools diagnostics in the QC log | Correct the input or explicitly configure the real field |
| Dataset sample cannot be selected | Dataset ID is not a VCF sample and the VCF has multiple samples | Run `bcftools query -l` | Supply the intended dataset/sample ID or a single-sample VCF |
| Frequency discordance is high | Allele orientation, build, ancestry, or reference release differs | Compare VCF metadata and reference provenance | Use a scientifically compatible reference; do not relax the threshold blindly |
| Resume refuses an output | A recorded output changed after validation | Read `checkpoint_events.log` and compare the recorded checksum | Preserve/review the file; use `--overwrite` only for intentional replacement |

## Limitations

QC summary assesses records but does not materialize a passing VCF, repair
alleles or coordinates, infer reference compatibility, or establish causality.
Independent rule counts overlap by design. The single extracted table reduces
repeated VCF scans but still requires disk space proportional to the queried
records and memory suitable for the streaming aggregations.

## Scientific references

- [bcftools query and statistics documentation](https://samtools.github.io/bcftools/bcftools.html)
- [Genome Reference Consortium GRCh37 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [GWAS-VCF specification](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
