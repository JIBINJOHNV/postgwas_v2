# Troubleshooting

Start with the first validation or execution error, not the final missing-file
symptom produced by a failed upstream stage.

## Command is not found

Confirm the environment is active and the checkout was installed. In a
container, run the command inside the built image. A working Python CLI does not
prove external tools such as bcftools, PLINK, MAGMA, GCTA, FINEMAP, or R are on
the executable path.

## Configuration is rejected

Validate the same file used for execution:

```console
postgwas config validate --config /absolute/path/to/run.yaml
```

Then inspect the relevant module's resolved configuration:

```console
postgwas config show \
  --config /absolute/path/to/run.yaml \
  --module MODULE
```

Do not add undocumented keys or copy defaults from `--help` when they disagree
with the exported YAML.

## Sample-sheet paths fail

Relative input and external-file paths are resolved from the sample sheet's
directory. Check the delimiter, header spelling, duplicate dataset IDs, empty
external files, and the mutually exclusive internal-versus-external frequency
and INFO mappings.

## Few variants remain after harmonisation or filtering

Inspect counts by rejection reason. Common causes include build mismatch,
coordinate parsing, invalid or swapped alleles, population-frequency mismatch,
palindromic-variant policy, INFO/MAF thresholds, sample-size failures, and
variant-ID joins. Do not loosen thresholds until the cause is understood.

## Reference joins are poor

Check build, chromosome naming, variant-ID representation, allele order,
population, resource release, compression, and indexes. Compare a small set of
variants manually before starting a full rerun.

## An external tool failed

Read both the PostGWAS canonical log and the external tool's own log. Record the
exact command and version. Confirm input files are non-empty and all required
companion files exist. Treat partially created outputs as incomplete.

## Pipeline order is surprising

The planner inserts registered dependencies and can repeat formatting around
imputation. Print or inspect the plan; do not infer order only from the target
list. Harmonisation and pathway enrichment are standalone-only.
