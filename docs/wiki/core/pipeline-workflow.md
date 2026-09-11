# Pipeline workflow

PostGWAS supports both standalone module commands and a pipeline command. Use a
standalone command when you already have that module's required input artifacts.
Use the pipeline when PostGWAS should run the registered preceding steps for one
or more final analyses. Pipeline mode is an orchestration choice, not a different
statistical method or an instruction to run every module.

For a complete connected example, start with [Quick Start](../getting-started/quick-start.md).
For the same MAGMA analysis with manually supplied intermediates, see the
[direct-mode comparison](running-modules-independently.md#example-magma-with-and-without-pipeline-orchestration).

## Harmonisation comes first

Raw summary statistics are prepared with the standalone harmonisation command:

```console
postgwas harmonisation --help
```

Harmonisation is not currently a pipeline-selectable step. Complete it first
when starting from raw data, then provide its validated study VCF to downstream
commands. An existing compatible PostGWAS-harmonised VCF can be used directly;
an arbitrary GWAS-VCF cannot. The entry check requires PostGWAS provenance,
the supported build declaration, one sample, an index, and the configured field
contract. Read [Input and Output Contracts](input-output-contracts.md).

Pathway enrichment is also standalone-only: it starts from a gene list, not the
pipeline study-VCF interface.

## Select final analyses

Pipeline targets describe the results you want, not every prerequisite you
expect to run. For example:

```console
postgwas pipeline \
  --modules finemap \
  --clumping-methods standard \
  --finemap-method susie \
  --help
```

PostGWAS validates the selected targets and shows their execution order. Shared
preceding steps normally run once. Formatting before imputation and formatting
the resulting VCF for downstream analysis are separate steps because they use
different data.

The contextual help page shows the steps and options for the selected methods.
Fine-mapping first asks for workflow choices when they are not explicit; include
both clumping and fine-mapping method choices to inspect that concrete plan.
Help inspects the interface and plan, not the contents of your study or references.

You supply the PostGWAS study VCF, dataset ID, output directory, external
references, and study-specific choices. The pipeline supplies intermediate
artifacts created by its selected upstream stages and hides those input options
from contextual help. It does not download missing reference panels or silently
replace them with a different build or ancestry.

## How the execution order is decided

The planner resolves method-dependent prerequisites, expands them recursively,
rejects unavailable targets or cycles, and produces a deterministic execution
order. There is no mandatory all-module chain. This table covers every public
pipeline target as a separate invocation without optional workflow switches.
Arrows give the scheduled order, not an exclusive input/output relationship
between adjacent stages.

For `finemap`, `caldera` and `flames`, the table explicitly selects
`--clumping-methods standard` and either `--finemap-method susie` or
`--finemap-method finemap`. The stage names are shared across these engines;
the formatter contracts, executables and reference requirements are not.

| Selected target and method | Execution order |
|---|---|
| `sumstat_filter` | `sumstat_filter` |
| `qc_summary` | `qc_summary` |
| `formatter` | `formatter` |
| `imputation` | `formatter → imputation` |
| `annot_ldblock` | `annot_ldblock` |
| `ld_clump`, `standard` only | `ld_clump` |
| `ld_clump`, `region` | `annot_ldblock → ld_clump` |
| `ld_clump`, `cojo-slct` | `formatter → ld_clump` |
| `finemap`, standard clumping | `ld_clump → formatter → finemap` |
| `magma` | `formatter → magma` |
| `gcta_gene` | `formatter → gcta_gene` |
| `gcta_cojo` | `formatter → gcta_cojo` |
| `magmacovar` | `formatter → magma → magmacovar` |
| `single_cell`, `magma_celltype` or `scdrs` | `formatter → magma → single_cell` |
| `single_cell`, `ldsc_celltype` only | `formatter → single_cell` |
| `pops` | `formatter → magma → pops` |
| `kpops` | `formatter → magma → kpops` |
| `caldera`, standard clumping | `formatter → ld_clump → magma → pops → finemap → caldera` |
| `flames`, standard clumping | `ld_clump → formatter → magma → magmacovar → pops → finemap → flames` |
| `heritability` | `formatter → heritability` |
| `mixer` | `formatter → mixer` |
| `manhattan` | `manhattan` |

FLAMES receives evidence from fine-mapping, MAGMA, MAGMAcovar and PoPS.
CALDERA consumes PoPS and credible sets, not K-POPS or MAGMAcovar outputs.
Its pipeline currently requires GRCh37; GRCh38 direct-mode support does not
imply GRCh38 pipeline support.
The planner currently schedules formatter before clumping for CALDERA and
after standard clumping for FLAMES because their prerequisite traversal differs.

The packaged clumping choice is `region` plus `standard`, so using it also
adds `annot_ldblock`; the standard-only rows above are explicit alternatives.
Harmonisation and pathway enrichment are not pipeline targets.
For multiple targets, inspect the combined plan rather than concatenating rows.

Combined clumping methods combine their prerequisites. Fine-mapping requires
standard clumping for its pipeline locus artifact and an engine-specific
formatter input. Selecting region pruning adds LD-block annotation; selecting
COJO clumping changes where its formatter input is needed. Do not use a fixed
diagram to infer the order of every fine-mapping combination.

Filtering, imputation, plotting, and QC appear only when the selected targets or
switches request them. The planner can repeat formatting around imputation
because the post-imputation data need new downstream inputs. Shared prerequisites
are otherwise scheduled once where compatible.

Imputation also re-harmonises its generated statistics internally, so selecting
it requires the full harmonisation resources as well as the PRED-LD reference.
The fact that public `harmonisation` is not a selectable pipeline target does
not remove those internal requirements.

Step directories are numbered by position in this resolved plan, not by a
permanent module number. A different target or method selection can therefore
change both directory numbers and required references.

## Optional workflow stages

The pipeline also accepts workflow switches that add filtering, imputation,
plotting, or heritability estimation without naming them as targets:
`--apply-filter`, `--apply-imputation`, `--apply-manhattan`, and
`--heritability`.

```console
postgwas pipeline \
  --modules magma \
  --apply-imputation \
  --apply-manhattan \
  --heritability \
  --help
```

When filtering and imputation are combined, the planner schedules filtering
before imputation and a distinct post-imputation filtering stage. Always inspect
the displayed plan rather than assuming that a command-line option maps to only
one physical step.

There is a current CLI-composition limitation: adding `--apply-filter` to
`ld_clump`, or to `finemap` including explicit standard/SuSiE selection, can
fail before analysis with an argparse conflict for `--remove-mhc` /
`--no-remove-mhc`. MAGMA with `--apply-filter` also fails during help assembly,
with a duplicate `--mhc-chrom` option. This is not a general prohibition on
combining filtering with any module that has a genome-build setting; the `gcta_cojo` and `qc_summary`
help combinations do construct successfully. For the affected routes, run
`postgwas sumstat_filter` separately and pass its validated filtered VCF to the
pipeline. This documentation workaround does not fix the parser defect.

## Export the matching configuration

Generate a configuration for a target and all required preceding steps:

```console
postgwas config export \
  --pipeline finemap \
  --style full \
  --output finemap_pipeline.yaml
```

Review the exported method selections and complete all required resource paths
before running the analysis. Continue to name the requested targets with
`--modules`; exporting settings is not a run, a resource download, or proof of
input readiness.

## Common input and resource validation

Before scientific stages run, the pipeline checks the entry GWAS-VCF and runs
the registered preflight for each selected module. It collects independent
module failures and presents the observed file checks in a common style.
The YAML audit is written at `pipeline.validation.report_file`; terminal detail
uses one section per file, with indexes grouped. The shared
reporter combines module-specific file facts into that section and shows only
new findings at later stages. Cross-file comparisons refer to stable file
numbers, while stage progress and scientific results remain separate. The
version-2 audit stores each exact file path once with distinct check observations.
The shared
`logging.file_validation.max_screen_files` defaults to 20 successful file
sections and never hides problems; additional successful files are summarized
on screen and retained in full in the saved record set. The audit is updated before
successful completion is announced
or when execution fails, with checks observed at later stages and the failure.

Existing checks have different coverage: some inspect complete file contents,
while others establish availability or header/index validity only. Files produced
by earlier pipeline stages are checked when their consumers can inspect them.
Read [Pipeline Input Validation](pipeline-input-validation.md) before treating
a startup pass as evidence of scientific compatibility.

## Execute the plan

A configured fine-mapping run has the following shape after the exported YAML
has been edited with all required references. This example explicitly selects
standard clumping and SuSiE, so the configuration must supply the prepared
pairwise-LD reference for standard clumping and a matching PLINK reference for
SuSiE, together with the compatible build/population settings:

```console
postgwas pipeline \
  --modules finemap \
  --clumping-methods standard \
  --finemap-method susie \
  --vcf study.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config finemap_pipeline.yaml
```

The pipeline stops at the failed step if validation, a scientific module, or an
external tool fails. A later step is not reported as complete unless execution
reaches it successfully.

## Related pages

- [Configuration](configuration.md)
- [Input and Output Contracts](input-output-contracts.md)
- [Pipeline Input Validation](pipeline-input-validation.md)
- [Logging and Reproducibility](logging-and-reproducibility.md)
