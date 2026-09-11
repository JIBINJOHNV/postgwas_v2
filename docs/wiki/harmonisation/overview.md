# Harmonisation

## Purpose

Harmonisation turns one or more study-specific GWAS summary-statistics tables
into allele-aware, build-resolved GWAS-VCF files for PostGWAS analyses.

## What the analysis does

For each sample-sheet row, PostGWAS validates the declared columns and files,
normalizes statistics and alleles, resolves genome build, processes chromosomes,
creates GWAS-VCF, annotates variants against configured references, creates the
opposite-build VCF by liftover, merges chromosome outputs, and writes QC,
rejection, provenance, and optional input-to-VCF concordance reports.

Study EAF requires exactly one declared source: a study column or an explicit
external file/column pair. INFO has a different rule: internal INFO takes
priority over an external source; when neither is supplied, the public command
requires an explicit `--fixed-info VALUE` choice. External proxy and fixed INFO
are recorded as such, not presented as study-measured imputation quality. The
comparison-frequency panel is used for validation and orientation checks; it
is not a silent replacement for missing study values.

## When to use it

Run harmonisation once for every raw GWAS dataset before using the downstream
pipeline. Re-run it when the raw dataset, column mapping, genome reference,
population frequency source, or scientifically material policy changes.

## Input requirements

- A version-2 sample sheet with one row per dataset.
- Every referenced summary-statistics or external data file.
- A resource directory containing the configured build-check, FASTA, dbSNP,
  annotation, chain, and allele-frequency resources.
- `bcftools`, `tabix`, and the configured bcftools `liftover` plugin.
- An optional run YAML for overrides; defaults come from the canonical
  harmonisation YAML.

See [Harmonisation Sample Sheet](sample-sheet.md) and
[Reference Resources](../reference/reference-resources.md).

## Direct mode

Run harmonisation directly with a reviewed sample sheet and complete reference
tree. Replace the paths below with your own files.

### Minimal example

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --resource-directory /absolute/path/to/resources \
  --output-directory /absolute/path/to/results
```

### Full example

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --dataset-id STUDY \
  --run-config harmonisation.yaml \
  --resource-directory /absolute/path/to/resources \
  --output-directory /absolute/path/to/results \
  --comparison-af-source ALFA \
  --comparison-af-column EUR \
  --threads 8 \
  --memory-gb 32 \
  --seed 10 \
  --keep_gwas2vcf_intermediate \
  --validate
```

Replace the comparison population, build resources, paths, and compute limits
with values appropriate to the dataset. `--validate` is optional and is off
unless requested in CLI or YAML. `--keep_gwas2vcf_intermediate` is also
optional; without it, a successful run deletes the raw merged adapter VCF after
final validation while retaining both annotated build-specific VCFs.

## Pipeline mode

Harmonisation is standalone-only: `harmonisation` is not a supported
`postgwas pipeline --modules` target. Run it first, inspect its reports and
validated final GWAS-VCF, then pass that VCF to a supported downstream pipeline.
The [connected Quick Start](../getting-started/quick-start.md) shows the complete
raw-data → harmonisation → QC → MAGMA route, including the required references.
Pipeline orchestration does not generate the sample sheet or download the
harmonisation reference tree for you.

## Parameters

`--sample-sheet` is required. `--dataset-id` selects one row; otherwise all
rows are processed. `--run-config` supplies YAML overrides.
`--resource-directory` is required unless `resources.root` is set in YAML.
`--output-directory`, comparison AF source/column, compute limits, seed, and
concordance validation are overrideable.

Export the current defaults instead of copying them from this page:

```console
postgwas config export \
  --module harmonisation \
  --style full \
  --output harmonisation.yaml
```

CLI values override user YAML; user YAML overrides packaged YAML. See
[Harmonisation Configuration](../../modules/harmonisation/configuration.md).

## Processing steps

The command and input declaration are preflighted before full reading. After
the data reveal the genome build and observed chromosomes, a second exact
resource preflight must pass before partition files or chromosome workers are
created. The dataset is then harmonised per chromosome, converted and annotated
as VCF, and finally merged and assessed.

**[How Harmonisation Processes Your Data](processing-order.md)** is the
step-by-step walkthrough of that sequence. It lists all seven stages and the 29
numbered steps in the exact order they run, using the same labels that appear in
the log files, and explains what each step reads, decides, and writes. Read it
before interpreting a harmonisation log or QC report.

## Outputs

Each dataset receives merged VCFs in both configured builds, indexes,
intermediate evidence, rejection reports, QC reports, a run manifest, and
dataset/chromosome logs. The raw same-build GWAS-to-VCF merge is validated as a
required intermediate, then deleted after successful finalization unless
`--keep_gwas2vcf_intermediate` (or the corresponding YAML policy) is enabled.
See
[Harmonisation Outputs and QC](outputs-and-qc.md).

The source-build and lifted-build VCFs can have different record counts. Under
the packaged `policies.vcf.liftover_swap: exclude` policy, successfully lifted
REF/ALT-swapped records are excluded separately from plugin lift failures.
Inspect both categories, along with earlier rejected variants, before choosing
the build-specific VCF for downstream analysis.

## QC and logs

Review dataset status, inferred build, chromosome failures, counts at each
transformation, rejected reasons, frequency concordance, liftover loss, VCF
metrics, virtual filter counts, effective sample-size summaries, and optional
input-to-VCF concordance. Post-merge QC is report-only; earlier harmonisation
steps can already have rejected variants, including unresolvable palindromic
variants. The merged VCF is not replaced by the virtual QC-passed subset. A
record can fail multiple QC rules, so their individual counts are not additive.

## Interpretation

A harmonised VCF is analysis-ready only for references and modules that match
its build, allele convention, population, and sample-size definition. Retain
the manifest and QC evidence with every downstream result.

## Common problems

- Incorrect column names or delimiters.
- Supplying both internal and external EAF, or supplying neither EAF source.
- Supplying no INFO source without explicitly choosing `--fixed-info`; when
  both internal and external INFO are listed, internal INFO takes priority and
  the ignored external source is reported.
- Case-control rows missing cases or controls.
- Build-check or chromosome resource files missing.
- Low build-match evidence or excessive liftover loss.
- Palindromic or allele-frequency-discordant variants.
- Missing external executables or the liftover plugin.

## Limitations

Automatic inference is evidence-based but cannot repair mislabeled study
metadata. Liftover cannot guarantee identical variant representation across
builds. Optional concordance reports unmatched representations rather than
using their count as a value-integrity failure; same-position checks cannot
prove two differently represented variants are biologically identical and are
therefore restricted to allele-unmatched, one-to-one, same-type records. The
output QC subset is virtual, not a second filtered VCF.

## Scientific references

- [GWAS-VCF specification and implementation](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
- [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
- [GWAS Catalog harmonisation methods](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics)
- [Mallard et al., effective sample size](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22609)
