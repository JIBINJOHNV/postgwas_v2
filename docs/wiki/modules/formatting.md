# Formatting

## Purpose

Formatting reads one harmonised GWAS-VCF once and creates validated,
tool-specific inputs for downstream PostGWAS modules. It can also write one
additional user-named table directly from CLI column options.

## What the analysis does

The formatter extracts a typed canonical table with bcftools, selects configured
variant IDs independently for each target, infers quantitative versus binary design from case/control count values,
applies each target's required-field and numeric checks, transforms statistics,
and writes outputs atomically. It can produce MAGMA, GCTA gene, SuSiE, FINEMAP,
PRED-LD, LDSC, and MiXeR inputs and an optional custom table in one run.

## When to use it

Use it when a downstream analysis requires a tabular contract rather than
GWAS-VCF. Pipeline planning inserts it before dependent analyses and may repeat
it after imputation.

## Input requirements

A harmonised, biallelic GWAS-VCF; dataset ID; output directory; `bcftools`; and
at least one built-in format from CLI/YAML or `--custom-output` with `--id`.
Required VCF fields depend on the chosen target and requested custom columns.

## Command

```console
postgwas formatter --vcf PATH \
  [--format FORMAT [FORMAT ...]] [--custom-output FILE] [options]
```

## Minimal example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma ldsc
```

### Custom CLI output example

The custom table is additional to any built-in formats and requires no user
YAML changes. The output columns follow the order in which their options appear
on the command line.

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma \
  --custom-output STUDY_for_my_tool.tsv \
  --id SNP \
  --chr CHR \
  --pos BP \
  --alt A1 \
  --ref A2 \
  --beta BETA \
  --se SE \
  --p P \
  --eaf EAF \
  --n N \
  --variant-id-type unique
```

## Full example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma gcta_gene susie finemap pred_ld ldsc mixer \
  --run-config formatting.yaml \
  --resume \
  --bcftools bcftools
```

## Parameters

`--format` accepts `magma`, `gcta_gene`, `susie`, `finemap`, `pred_ld`, `ldsc`,
and `mixer` plus documented aliases. When omitted, formats come from YAML.
`--variant-id-type rsid` extracts an rsID from the VCF `ID` field;
`--variant-id-type unique` constructs the configured chromosome-position-REF-ALT
identifier. `variant_identifiers.target_types` can set different conventions
for different outputs. MAGMA pipeline mode inspects BIM field 2 and sets only
the MAGMA target automatically.

`--custom-output FILE` activates the additional custom table. `--id NAME` is
mandatory and uses the same `--variant-id-type rsid|unique` policy. Each other
option requests one field and supplies its output header:

| Option | Value written |
|---|---|
| `--id NAME` | Selected rsID or configured chromosome-position-REF-ALT unique ID |
| `--chr NAME` | Normalized chromosome |
| `--pos NAME` | One-based position |
| `--ref NAME` | Other/non-effect allele |
| `--alt NAME` | Effect allele |
| `--beta NAME` | Effect estimate |
| `--se NAME` | Standard error |
| `--z NAME` | Z statistic |
| `--lp NAME` | Input `-log10(P)` |
| `--p NAME` | Raw P converted from `LP`, bounded by the configured minimum |
| `--eaf NAME` | Effect-allele frequency |
| `--maf NAME` | `min(EAF, 1-EAF)` |
| `--n NAME` | Total sample size |
| `--neff NAME` | Effective sample size |
| `--n-case NAME` | Case sample size |
| `--n-control NAME` | Control sample size |
| `--info NAME` | Imputation INFO score |

There is no permissive missing-value option. Every requested field is required
for every written row; rows with missing, non-finite, or scientifically invalid
requested values are excluded and counted. Output header names must be unique.

Resume reuses outputs only after provenance, content, configuration, and
freshness validation. If a completed run directory is copied, resume validates
the copied files at the currently configured destinations, rebases every
returned artifact and manifest path to that directory, and rejects a recorded
artifact whose relative filename does not match the current configuration.
`--overwrite` takes precedence. Output schemas,
filenames, chromosomes, minimum representable P value, and MiXeR QC come from
the generated [Configuration Defaults](../reference/configuration-defaults.md).

## Processing steps

The formatter validates selection, renders every selected output destination,
and rejects filename collisions before VCF extraction. This preflight includes
named outputs and every configured chromosome partition. It then optionally
validates a completion manifest, extracts `N_ALT=1` records once, types numeric
fields, selects each target's configured ID convention, and runs each exporter
in canonical order. The custom exporter reuses the same one-pass extraction and
does not trigger study-design inference. Study design is inferred only when
LDSC or MiXeR is selected,
because only those formatter contracts interpret sample size differently for
binary and quantitative traits. MAGMA, SuSiE, and FINEMAP validate only their
own configured fields and do not require case/control columns. Rows missing or
violating a target's required fields are excluded for that target and counted.
Before PRED-LD field validation, rows whose normalized chromosome is absent
from the configured `chromosomes` list are excluded and counted by chromosome
label.

## Outputs

- MAGMA: `<dataset>_magma_snp_loc.tsv` and `<dataset>_magma_p_values.tsv`.
- GCTA gene: `<dataset>_gcta.ma`, shared by fastBAT and mBAT-combo.
- SuSiE: `<dataset>_susie.tsv`.
- FINEMAP: `<dataset>_finemap.tsv`.
- PRED-LD: `pred_ld/<dataset>_chr<chromosome>_pred_ld_input.tsv` for configured chromosomes.
- LDSC: `<dataset>_ldsc_input.tsv`.
- MiXeR: `<dataset>_mixer.sumstats.gz`.
- Custom: the relative filename supplied to `--custom-output`.
- Provenance: `logs/<dataset>_formatter.log`,
  `run_metadata/resolved_config.yaml`, and
  `run_metadata/formatter_completion.yaml`.

## QC and logs

For each target, review the selected ID type, identifier exclusions, rows
in/out, exclusions for missing or invalid values,
P values bounded at the configured numeric minimum, written schema, inferred
study design when required, sample-size mode, and output fingerprints. PRED-LD additionally
records `rows_excluded_unconfigured_chromosome` in its result and lists each
excluded normalized chromosome with its row count in the canonical log.
The custom result additionally records its ordered field roles, output headers,
identifier convention, and missing/invalid requested-field exclusions.

## Interpretation

Allele direction is preserved: GWAS-VCF ALT is the effect allele. FINEMAP gets
minor allele frequency; MAGMA and GCTA mBAT receive total sample size; SuSiE
uses effective N; LDSC uses case/control counts for binary traits or N from NCO
for quantitative traits; MiXeR receives effective N.

## Common problems

No format selected, colliding resolved output filenames, no usable IDs of the
selected type, duplicate selected IDs, missing VCF fields, incomplete
per-variant case/control counts, stale outputs, a changed resolved config,
manifest artifact paths that do not match the current configured filenames, or
zero rows satisfying one target's scientific contract.

## Limitations

Formatting does not run the external analyses and cannot make an inappropriate
reference panel compatible. A custom table guarantees the requested formatter
semantics but does not assert that the table satisfies an arbitrary external
tool's complete input contract. When LDSC or MiXeR requires study design, it is
inferred from count values rather than a separate VCF trait declaration.

## Scientific references

See [Scientific References](../reference/scientific-references.md) and the
documentation for the downstream tool receiving each export.
