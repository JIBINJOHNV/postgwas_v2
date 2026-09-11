# Imputation

## Purpose

The current imputation command uses PRED-LD to infer missing summary statistics
from formatted per-chromosome inputs and an LD reference, then post-processes
and re-harmonises the available results.

## What the analysis does

For each configured chromosome, PRED-LD selects LD-linked reference variants,
propagates effect direction using LD correlation, and produces imputed and
observed records. PostGWAS combines chromosomes, reports correlations for
duplicated imputed/observed markers, propagates sample metadata from LD-linked
donors, writes a version-2 harmonisation sample sheet, and invokes
harmonisation on the combined table.

## When to use it

Use it only when missing summary statistics materially limit a downstream
analysis and a build- and ancestry-matched PRED-LD resource exists. Compare
imputed and observed statistics before relying on the expanded dataset.

## Input requirements

- Formatter-created `pred_ld/` files named
  `<dataset>_chr<chromosome>_pred_ld_input.tsv`.
- A compatible PRED-LD LD-reference directory and PostGWAS resource directory.
- An explicit genome-build declaration for the PRED-LD release. PRED-LD resource
  filenames do not encode a build, so PostGWAS cannot safely infer it.
- Population, reference mode, dataset ID, output directory, and adequate memory.
- Full harmonisation prerequisites because imputed output is re-harmonised.

Pipeline preflight validates the harmonised GWAS-VCF once, compares its build
with `modules.imputation.genome_build`, validates every chromosome-specific
`TOP_LD` reference file and the harmonisation software, and only then allows the
formatter to create PRED-LD input tables. The validated resource identities are
rechecked before imputation so a replaced reference cannot be consumed without
restarting validation.

## Command

```console
postgwas imputation --pred-ld-input-directory PATH [options]
```

## Minimal example

```console
postgwas imputation \
  --pred-ld-input-directory formatted/pred_ld \
  --imputation-ld-reference reference/pred_ld \
  --resource-directory resources \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas imputation \
  --pred-ld-input-directory formatted/pred_ld \
  --imputation-engine pred_ld \
  --imputation-ld-reference reference/pred_ld \
  --imputation-r2-threshold 0.8 \
  --imputation-minimum-maf 0.001 \
  --ref TOP_LD \
  --resource-directory resources \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results \
  --threads 2
```

## Parameters

Only `pred_ld` is implemented in direct mode. PostGWAS currently supports its
`TOP_LD` handoff because the bundled post-processor consumes the corresponding
`LD_info_TOP_LD` output contract. Canonical `modules.imputation` YAML supplies
minimum LD r² 0.8, reference MAF 0.001, population EUR, reference mode
`TOP_LD`, and Pearson correlation. The reference directory and genome build
are intentionally unset and required from YAML or CLI. CLI options override
matching YAML values. The runner applies additional memory-aware scheduling;
begin with conservative parallel settings for a new system.

Worker concurrency is bounded by the resolved `execution.memory_gb` budget
(`--memory-gb` overrides it), not host RAM. Resource preflight stops before
formatter or PRED-LD execution when that budget cannot fit the configured
`memory_gb_per_worker` reservation. The runner repeats this guard before
creating its outputs and logs both the budget and selected worker count.
Reservations are planning estimates, not enforced operating-system limits.

`--correlation-method` controls the observed-versus-imputed post-processing
correlation and overrides `engines.pred_ld.correlation_method`. Worker memory,
free-memory thresholds, polling intervals, large-chromosome classification, and
the preferred chromosome schedule are controlled under
`modules.imputation.engines.pred_ld` in YAML. These scheduling settings affect
resource use and execution order, not the statistical model. The same section's
`reference_subdirectory_template`, `reference_file_template`, and
`reference_file_kinds` describe the installed PRED-LD release layout, so the
resource path convention is validated configuration rather than Python logic.

## Processing steps

Read the prevalidated LD reference in place, schedule chromosome workers, run
the bundled PRED-LD entry point, require both PRED-LD outputs from every
configured chromosome, combine logs, join imputed variants to donor statistics,
compute duplicate beta and Z correlations using the configured method, create
combined and correlation tables, archive raw chromosome files, generate the
harmonisation sheet, and re-harmonise to a sibling `imputed_harmonised`
directory. Reading the reference in place avoids copying the complete panel for
each run.

## Outputs

Intermediate and post-processing outputs include:

- `imputation_results_chr<chromosome>.txt` and
  `LD_info_TOP_LD_chr<chromosome>.txt` before consolidation;
- `<dataset>_master.log` and `<dataset>_predld_combined.log`;
- `<dataset>_PREDLD_allchr.tsv.gz`;
- `<dataset>_PREDLD_correlations.tsv`;
- `<dataset>_harmonisation_sample_sheet.csv`;
- `<dataset>_imputation_results.txt.gz` and
  `<dataset>_LD_info_TOP_LD.txt.gz` archives;
- re-harmonised GWAS-VCFs and QC under `imputed_harmonised`.

## QC and logs

Verify every expected chromosome succeeded, inspect imputed-marker counts and
duplicate beta/Z correlations, confirm donor-derived NC/SS/AF/SI completeness,
and review the entire second harmonisation QC. A missing, empty, or failed
configured chromosome now fails the PRED-LD stage; the combined log records all
affected chromosomes and their reasons.

## Interpretation

Imputed statistics are model- and reference-dependent estimates, not observed
study results. LD population, build, allele orientation, donor coverage, and r²
threshold directly affect their validity.

## Common problems

Missing chromosome formatter files, incompatible reference build/population,
insufficient memory, missing PRED-LD side files, no successful chromosome,
incomplete donor statistics, or failure of the re-harmonisation prerequisites.

## Limitations

Only PRED-LD `TOP_LD` is currently supported. The supplied reference panel's
genome build must be declared explicitly because its filenames do not encode
that provenance.

## Scientific references

- [PRED-LD software and method documentation](https://github.com/pbagos/PRED-LD)
