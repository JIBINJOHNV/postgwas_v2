# Fine Mapping

## Purpose

Fine mapping estimates which variants within selected association loci are
most compatible with causal effects using SuSiE-RSS or FINEMAP and matched LD.

## What the analysis does

PostGWAS constructs range- or point-based locus windows, applies an LP inclusion
threshold and optional MHC exclusion, aligns formatted summary statistics with
a PLINK reference, computes locus LD, runs the selected model in parallel, and
combines successful locus results into QC and FLAMES-compatible outputs.

## When to use it

Use it after association loci have been defined and after producing the
method-specific formatter table. Fine mapping is meaningful only when LD,
build, ancestry, alleles, sample size, and locus variants are compatible.

## Input requirements

- A locus TSV: `CHR START END LP` for `range`, or `CHR POS LP` for `point`.
- `<dataset>_susie.tsv` for SuSiE or `<dataset>_finemap.tsv` for FINEMAP.
- A matching PLINK BED/BIM/FAM reference prefix.
- SuSiE: PLINK plus the required R/susieR environment.
- FINEMAP: PLINK 2, bgenix, LDstore, and FINEMAP executables.

## Command

```console
postgwas finemap --finemap-method {susie,finemap} [options]
```

Pipeline help is method-aware. First display the available selectors:

```console
postgwas pipeline --modules finemap --help
```

Then request the complete options for one valid workflow, for example:

```console
postgwas pipeline --modules finemap --clumping-methods standard --finemap-method susie --help
```

Fine-mapping currently requires `standard` among the selected clumping methods
because its genomic-risk-locus table supplies the model boundaries. `region`
and `cojo-slct` may be added for parallel locus summaries but do not replace
that locus source. With `standard` alone, the pipeline runs standard LD
clumping, creates the engine-specific formatter table, and then fine-maps the
resulting loci. LD-block annotation and its BED resources are added only when
`region` is selected. When `cojo-slct` is selected, formatter runs before LD
clumping because the COJO analysis consumes its GCTA table.

## Locus boundaries

The canonical boundary flank is 500 kb on each side. In `point` mode, a locus
at `POS` therefore spans `POS - 500 kb` through `POS + 500 kb` by default. In
`range` mode, PostGWAS expands the supplied or pipeline-generated interval from
`START - 500 kb` through `END + 500 kb`. This is a configurable operational
default, not a universal biological boundary; the LD architecture and study
design still determine whether a different region is justified.

Direct mode accepts predefined loci through `--locus-file`. Use
`--locus-type range --window-kb 0` when `START` and `END` are already the exact
model boundaries, or override the flank with `--window-kb INT`. Pipeline mode
does not accept a separate locus-file override: it uses the validated
genomic-risk-locus table produced by its required standard-clumping step, then
applies the resolved `--window-kb`/`locus_window_kb` flank.

## Minimal example

```console
postgwas finemap \
  --finemap-method susie \
  --susie-input-file formatted/STUDY_susie.tsv \
  --locus-file loci.tsv \
  --locus-type range \
  --window-kb 0 \
  --finemap-ld-reference reference/1000G_EUR \
  --plink plink \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas finemap \
  --finemap-method finemap \
  --finemap-in-files formatted/STUDY_finemap.tsv \
  --locus-file loci.tsv \
  --locus-type range \
  --window-kb 500 \
  --lp-threshold 7.3 \
  --finemap-ld-reference reference/1000G_EUR \
  --plink plink2 \
  --n-causal-snps 5 \
  --prob-cred-set 0.95 \
  --sss \
  --minimum-memory-per-worker-gb 14 \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

The canonical YAML defaults to SuSiE, range loci, 500-kb flanks, LP 7.3, 14 GB per
worker, and skipping chr6:25–35 Mb. SuSiE defaults to 10 components and 180-s
LD/model timeouts. FINEMAP defaults to SSS, 100,000 iterations, five causal
SNPs, 0.95 credible-set probability, and the other values shown by live help.
CLI help displays these YAML values, and only explicitly supplied CLI options
override them.

## Processing steps

First validate all summary-statistic rows and locus definitions. Then validate
the PLINK BED/BIM/FAM dimensions, BED binary consistency, locus-variant BIM
ID/coordinate concordance, per-locus reference coverage, and every
engine-specific executable. The screen reports the measured variant, locus,
reference, overlap, and runtime evidence before locus preparation begins. The
workflow then harmonizes summary and reference alleles, computes and QCs LD
matrices, schedules loci within memory limits, fits the model with timeouts,
records each locus status, merges successful outputs, and creates minimal
FLAMES handoff files.

## Outputs

Both engines use the same categorized layout: published tables, primary
credible sets, and plots are under `results/`; status tables, warnings, and
scientific audits are under `quality_control/`; configuration, software
versions, and run logs are under `run_metadata/`; and the authoritative
post-overlap handoff is under `downstream_inputs/flames/`. Caches, fitted model
objects, prepared locus inputs, primary handoff files, and engine workspaces are
isolated under `intermediate_files/`. Successful SuSiE worker copies are
removed after the published aggregates validate.

Every completed analysis also writes
`results/<dataset>_fine_mapping_report.html`. This self-contained report starts
with the authoritative post-overlap result, but it also retains the complete
primary-locus outcome history, input/reference preflight evidence, overlap and
joint-rerun provenance, engine recovery or causal-count model-selection audit,
resolved scientific parameters, software versions, and links to the durable
source files. Its credible-set and variant tables are searchable, sortable,
and paginated in the browser; all rows remain embedded rather than being
reduced to a significance threshold. The report fails validation instead of
publishing if a final set file has invalid posterior values, inconsistent
coverage metadata, or a member count that disagrees with its manifest.

## QC and logs

Review attempted/successful/failed loci, variants matched and removed, allele
flips, LD asymmetry/diagonal/eigenvalue checks, model convergence, credible-set
purity/coverage, timeouts, and software versions. The run can complete with a
mixture of successful and failed loci; never infer complete locus coverage from
the top-level status alone.

`quality_control/input_and_resource_validation.tsv` preserves the run-level
input/reference/runtime checks displayed in the first two progress stages.

## Interpretation

Posterior inclusion probabilities and credible sets are conditional on the
chosen locus, variant set, LD matrix, prior/model, and maximum causal effects.
They do not prove biological causality.

For SuSiE, set membership is defined by one single-effect component's `alpha`
values at the configured coverage target, while the reported member value is
the model-wide PIP. Its sum inside a set is therefore not the component's
coverage. For FINEMAP, the member value is the SNP posterior probability in the
selected model; the model's `Post-Pr` is a separate probability used to select
the causal-count model. The HTML report labels these quantities separately.

## Common problems

Low summary/reference overlap, allele mismatch, wrong build/population, singular
or non-positive LD, insufficient memory, timeout, no locus passing LP, MHC-only
loci, or missing external executables.

## Limitations

Partial locus success requires manual review. Configuration is unified through
the canonical YAML; explicitly supplied CLI values override matching YAML
fields and are recorded in the resolved run configuration.

## Scientific references

- [Wang et al. 2020, SuSiE](https://doi.org/10.1111/rssb.12388)
- [SuSiE-RSS documentation](https://stephenslab.github.io/susieR/reference/susie_rss.html)
- [Benner et al. 2016, FINEMAP](https://doi.org/10.1093/bioinformatics/btw018)
- [FINEMAP documentation](https://christianbenner.com/)
- [Example SuSiE study using 500-kb lead-variant windows](https://pmc.ncbi.nlm.nih.gov/articles/PMC11293972/)
