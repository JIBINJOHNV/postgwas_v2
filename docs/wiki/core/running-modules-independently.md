# Running Modules Independently

Most analysis commands can be invoked directly, but direct execution does not
remove their upstream input requirements. Pipeline mode supplies outputs from
preceding stages, but it does not make incompatible inputs scientifically valid.

## Standalone preparation

`harmonisation` is intentionally standalone and precedes the downstream
pipeline. Use it to convert mapped raw summary-statistics columns into the
harmonised GWAS-VCF and associated QC evidence.

## Downstream standalone commands

Public analysis commands include filtering, formatting, imputation, LD
annotation, LD clumping, LDSC, fine-mapping, MAGMA, MAGMAcovar, GCTA gene
analysis, single-cell integration, PoPS, FLAMES, MiXeR, Manhattan plotting,
pathway enrichment, and QC. Their exact installed names are shown by
`postgwas --help`.

Before running one directly:

1. Read its `--help` output from the installed checkout.
2. Confirm the required upstream artifact exists and was produced with matching
   build, ancestry, alleles, and sample-size definitions.
3. Inspect the module's resolved YAML configuration.
4. Verify external executables and reference files.
5. Use a new output directory or dataset identifier to avoid confusing outputs
   from different runs.
6. Retain the command, configuration, logs, and QC summaries.

## Important module relationships

Important relationships include formatting before imputation, MAGMA, GCTA gene analysis, LDSC,
fine-mapping, and MiXeR; LD annotation before LD clumping; LD clumping plus
formatting before fine-mapping; MAGMA before MAGMAcovar and PoPS; and
fine-mapping, MAGMAcovar, and PoPS before FLAMES.

The planner can repeat formatting around imputation because the imputed data
must be exported again for downstream consumers. Inspect the printed plan rather
than assuming each module runs once.

MAGMA cell typing consumes an exact MAGMA `.genes.raw` file and an aggregated
single-cell expression covariate matrix in direct mode. Pipeline mode supplies
the gene results through the MAGMA dependency, but the biological covariate
matrix remains an explicit user input.

Pathway enrichment and harmonisation are standalone-only.
