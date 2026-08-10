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
- Population, reference mode, dataset ID, output directory, and adequate memory.
- Full harmonisation prerequisites because imputed output is re-harmonised.

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
  --correlation-method pearson \
  --resource-directory resources \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results \
  --threads 2
```

## Parameters

Only `pred_ld` is implemented in direct mode. The direct defaults are minimum
LD r² 0.8, reference MAF 0.001, population EUR, reference mode `TOP_LD`, and
Pearson post-processing correlation. The runner applies additional memory-
aware scheduling; begin with conservative parallel settings for a new system.

## Processing steps

Copy the LD reference into a temporary workspace, schedule chromosome workers,
run the bundled PRED-LD entry point, validate chromosome output presence,
combine logs, join imputed variants to donor statistics, compute duplicate beta
and Z correlations, create combined and correlation tables, archive raw
chromosome files, generate the harmonisation sheet, and re-harmonise to a
sibling `imputed_harmonised` directory.

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
and review the entire second harmonisation QC. The combined log is essential
because the worker return value can be true when only a subset succeeds.

## Interpretation

Imputed statistics are model- and reference-dependent estimates, not observed
study results. LD population, build, allele orientation, donor coverage, and r²
threshold directly affect their validity.

## Common problems

Missing chromosome formatter files, incompatible reference build/population,
insufficient memory, missing PRED-LD side files, no successful chromosome,
incomplete donor statistics, or failure of the re-harmonisation prerequisites.

## Limitations

Alternative engines listed in YAML are not implemented by the public direct
command. Current partial-chromosome handling requires manual coverage review,
and some scheduling/resource details are not canonical YAML settings.

## Scientific references

- [PRED-LD software and method documentation](https://github.com/pbagos/PRED-LD)
