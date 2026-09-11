# Migrating from v1

This guide translates examples in the [legacy PostGWAS README](https://github.com/JIBINJOHNV/postgwas#readme)
to the current `postgwas_v2` repository. Here, “v1” identifies that older
repository and its examples, not the `config_version` of a sample sheet.
There is no automatic migration step described here, and matching command
names do not imply identical defaults, retained variants or analysis results.

## 1. Install the current checkout separately

Follow [Installation](installation.md) for a fresh Mamba environment or a local
Linux-amd64 Docker build. Keep the old environment, inputs, commands and results
until you have compared the new run. The legacy `jibinjv/postgwas:1.3` image is
not the current v2 build; installing a new checkout does not update an existing
container image. The current Docker instructions use a locally built image,
not a replacement public image tag.

Reference files remain separate from software installation. Check their
build, population, allele/variant identifiers and gene identifiers against the
selected v2 method; an old resource directory is not automatically compatible.
See [Resource Setup](resource-setup.md).

## 2. Update scripts using the current public options

Common command names remain, including `harmonisation`, `formatter`, `magma`
and `pipeline`. Several option names in the legacy examples changed:

| Legacy example | Current interface | Scope |
|---|---|---|
| `--nthreads` | `--threads` | Commands exposing the shared compute settings |
| `--max-mem 16G` | `--memory-gb 16` | Numeric GB budget; do not append `G` |
| `--sample_id` | `--dataset-id` | Dataset selection or output identity, as documented by the command |
| `--outdir` | `--output-directory` | Output root |
| `--config` for the harmonisation CSV | `--sample-sheet` | Harmonisation study mappings |
| `--defaults` for the harmonisation YAML | `--run-config` | Harmonisation analysis/resource settings |
| `--snp_loc_file`, `--pval_file` | `--snp-location-file`, `--p-value-file` | Direct MAGMA inputs |
| `--ld_ref`, `--gene_loc_file` | `--magma-ld-reference`, `--gene-location-file` | The MAGMA example, not a global LD-reference alias |

Inspect help from the environment you will run:

```console
postgwas harmonisation --help
postgwas magma --help
postgwas pipeline --modules magma --help
```

Do not mechanically replace underscores in every old option. Inputs and
prerequisites depend on the selected method. For example, direct QC is `qc`,
whereas its pipeline target is `qc_summary`. Harmonisation and pathway
enrichment are standalone-only. Read [Pipeline Workflow](../core/pipeline-workflow.md)
and the [direct-mode guide](../core/running-modules-independently.md) before
reusing an old multi-module command.

## 3. Recreate and review the sample sheet

Replace the old `create_sumstat_map_pl.py` invocation with the current draft
generator, using a directory containing your original summary-statistics files:

```console
python -m postgwas.modules.harmonisation.sample_sheet_generator \
  --input-directory /path/to/summary-statistics \
  --output studies_v2.csv
```

This reads headers to propose mappings; it does not convert an old sample
sheet or prove that a study is ready. Review missing fields, mapping warnings
and the companion rejected-files report. Never invent counts, frequencies or
quality measurements to complete the draft.

The current sheet requires `config_version` equal to `2` and rejects unknown
fields. Names include `input_file` instead of `sumstat_file`, `dataset_id`
instead of `gwas_outputname`, and `chromosome_column`/`position_column` instead
of `chr_col`/`pos_col`. Regenerate from the source headers or copy a maintained
template; simply adding a version field to the old CSV is insufficient.

The templates are `examples/configs/harmonisation/sample_sheet_quantitative.csv`
and `examples/configs/harmonisation/sample_sheet_case_control.csv` in the
checkout. Replace their placeholder paths, column names and study values.
Declare the actual trait, effect and P-value types when known. Quantitative
total N uses `control_count` or `control_count_column`; case-control studies
need both case and control sources.

Review these changes in particular:

- Choose exactly one internal or external study-frequency source. An external
  source requires both `external_eaf_file` and `external_eaf_column`.
- External per-chromosome EAF/INFO files need an explicit `{chromosome}` path
  template, optionally with `{build}`. Do not reuse a bare prefix that relied
  on the old automatic chromosome suffix.
- INFO uses the current source priority. With no internal or external source,
  an explicit `--fixed-info` choice is required and recorded as assigned,
  not measured, quality.
- `resource_folder`, `output_folder` and `liftover` are not v2 sheet fields.
  Set resource/output locations through current CLI/YAML settings and review
  the current build-conversion policies rather than copying an old Yes/No flag.

Use the [complete sample-sheet contract](../harmonisation/sample-sheet.md) for
the accepted fields, source rules and relative-path behavior.

## 4. Rebuild YAML settings separately

The sample sheet maps raw study columns; YAML controls analysis policies,
resources and execution. Export a current template instead of assuming an old
defaults file has compatible keys:

```console
postgwas config export --module harmonisation --output harmonisation_v2.yaml
postgwas config validate --config harmonisation_v2.yaml
```

Edit the export with your verified settings and resource paths, then validate
it again. Configuration validation is not a check of every reference file.
Harmonisation reads the sheet through `--sample-sheet` and this YAML through
`--run-config`; keep them as separate files. Most direct commands accept run
configuration, but direct Manhattan and pathway enrichment currently do not.
See [Configuration](../core/configuration.md) and follow the connected
[Quick Start](quick-start.md) for complete analysis commands.

## 5. Use a new output root and compare the results

Do not reuse the old result directory or copy old checkpoint manifests into a
new run. Current resume requires matching checkpoints, inputs, settings,
software and validated outputs; file existence alone is insufficient. Review
the [resume and overwrite policy](../core/configuration.md#resume-and-overwrite)
before restarting an existing v2 run.

Update scripts that assume the legacy `00_harmonised_sumstat` directory.
The packaged harmonisation layout is now `<dataset>/harmonisation/` below the
chosen output root; pipeline stage numbers depend on the resolved plan.
See [Output Structure](../reference/output-structure.md).

An old VCF can be reused only if it satisfies the selected consumer's current
contract. If required PostGWAS provenance is missing, regenerate from the
original summary statistics; do not add headers manually. See
[Input and Output Contracts](../core/input-output-contracts.md).
Compare variant retention and exclusions, allele orientation, statistic and
sample-size meanings, gene mappings, reference versions and resolved settings
before comparing final estimates. A successful v2 command does not establish
numerical equivalence to a v1 analysis.
