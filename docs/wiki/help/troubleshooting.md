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
external files, and the exact `{build}`/`{chromosome}` template spelling.
Study EAF requires exactly one internal or external source. INFO is not mutually
exclusive: internal INFO takes priority over external INFO, and absent INFO
requires an explicit `--fixed-info` choice. Do not invent measured quality to
satisfy preflight. See the [sample-sheet contract](../harmonisation/sample-sheet.md#eaf-and-info-mappings).

## Software is installed but reference files are missing

The installer and its tool checks do not install the full harmonisation resource
tree or all downstream panels. Follow [Resource Setup](../getting-started/resource-setup.md)
for upstream acquisition, preparation gaps and required companion files. A
directory with the expected name is not sufficient. Keep reference generation
separate from the analysis; do not edit a reference during a running job.

## Few variants remain after harmonisation or filtering

First identify where counts changed: harmonisation rejection, lift failure,
post-liftover exclusion, explicit filtering, or downstream reference matching.
Inspect the harmonisation `rejected/` records, reason matrix and HTML report;
then inspect the selected downstream module's exclusions. The packaged
liftover policy excludes successfully lifted REF/ALT-swapped records separately
from failed lifts. Post-merge QC and the `qc` command report rule failures but
do not themselves remove variants; INFO/MAF filtering requires the filtering
step. A variant may fail more than one QC rule, so do not sum overlapping counts.

Common causes of actual losses include build mismatch, coordinate parsing,
unresolvable allele orientation, incompatible frequency references,
palindromic-variant policy, sample-size failures, and variant-ID joins. Do not
loosen thresholds until the cause is understood. See
[harmonisation outputs and QC](../harmonisation/outputs-and-qc.md).

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
