# Quick Start

This walkthrough takes one study from raw summary statistics to a QC assessment
and positional MAGMA gene-association results. It uses the same dataset name and
output paths throughout; no example data or references are downloaded by these
analysis commands.

```text
Raw summary statistics + studies.csv
  → harmonisation → PostGWAS GWAS-VCF
                      ├─ qc → report-only assessment
                      └─ pipeline --modules magma → formatter → MAGMA
```

The example downstream analysis uses GRCh37, a European LD reference, and the
NCBI37.3 Entrez gene-location file, matching the packaged positional MAGMA
configuration. Use it only when these choices are appropriate for the study.
Other builds, ancestries, and gene annotations need their matching references
and configuration; do not relabel a reference to make it fit the example.

## 1. Install and prepare the inputs

Complete [Installation](installation.md), then activate the local environment
or enter the documented Docker shell with input/reference/output mounts. For
Docker, use paths inside the container in all commands and sample-sheet rows.
Before running this example, replace the following paths with your real inputs:

| Example path | What must be present |
|---|---|
| `studies.csv` | A reviewed version-2 sample sheet containing a row whose `dataset_id` is `STUDY` and whose input path points to the real study file |
| `reference/harmonisation/` | The configured build-check, FASTA, dbSNP, frequency, annotation, and chain resources for the input and output builds |
| `reference/GRCh37_EUR_reference` | A compatible PLINK prefix with `.bed`, `.bim`, and `.fam` companions; the BIM identifiers must match the formatter's selected identifier convention |
| `reference/NCBI37.3.gene.loc` | The NCBI37.3 MAGMA gene-location annotation with Entrez gene IDs |

The [resource setup guide](resource-setup.md) separates upstream downloads from
PostGWAS-specific preparation and gives a
[positional MAGMA acquisition route](resource-setup.md#prepare-the-positional-magma-example).
These are placeholders, not a bundled demonstration dataset. The complete
harmonisation resource tree still needs manual preparation; stop here if any
required component is unavailable.

## 2. Prepare the sample sheet

Use one row per study and map the actual source columns. Declare the coordinate
and allele columns, effect estimate or Z score, P-value representation, exactly
one study-frequency source, and the real sample-size sources. Quantitative and
case-control studies have different sample-size requirements. Do not invent
case counts, control counts, allele frequencies, or imputation-quality values.

Check the study publication or data dictionary before accepting `auto` values:
set `trait_type` to the documented study design, `effect_type` to `beta` or
`odds_ratio` when known, and `p_value_type` to `raw` or `neglog10` when known.
Confirm that the effect allele, effect estimate and frequency describe the same
allele, and whether SE is on the beta/log-odds scale or the raw odds-ratio scale.
Header recognition does not establish these meanings. A case-control study's
total N alone is not enough to reconstruct its effective N; supply the real
cases and controls or
review the [documented Neff-only limitation](../harmonisation/sample-sheet.md#a-case-control-file-containing-only-precomputed-neff).

Repository templates are available at:

```text
examples/configs/harmonisation/sample_sheet_quantitative.csv
examples/configs/harmonisation/sample_sheet_case_control.csv
```

You can instead [generate a draft from the study headers](../harmonisation/sample-sheet.md#generate-a-draft-from-gwas-headers)
and review all mappings and missing fields. Successful draft generation is not
proof of analysis readiness.

Copy and complete the appropriate template, or edit the generated `studies.csv`;
change the selected row's dataset ID to
`STUDY`, and replace all example values. Relative input paths are resolved from
the sample sheet's directory. Follow the [sample-sheet contract](../harmonisation/sample-sheet.md)
for INFO-source choices, external sources, missing statistics, and accepted
column combinations. The templates map study INFO; if your study has no internal
or external INFO source, the command requires an explicit `--fixed-info` choice,
recorded as user-assigned rather than measured quality. Do not leave an INFO
column mapped when it does not exist. A sample sheet is not a YAML run configuration.

## 3. Harmonise the study

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --dataset-id STUDY \
  --resource-directory reference/harmonisation \
  --output-directory results/harmonisation \
  --validate
```

Harmonisation is a standalone preparation command, not a pipeline target. It
validates and prepares each dataset, processes chromosomes, then merges and
assesses the results. It can reject variants; it does not promise to preserve
every input row. Build conversion can also lose records through failed liftover
or configured exclusions after a successful lift.
`--dataset-id STUDY` selects that sample-sheet row; omit it to process all rows.
The optional `--validate` used here compares the original study with its
same-build final VCF and records concordance evidence after harmonisation.

### Review the decisions used by this example

The commands use packaged policies unless you explicitly override them. The
following choices affect interpretation, not just file layout:

| Decision | Packaged behaviour to review |
|---|---|
| Frequency reference | ALFA EUR is used by both the default tabular comparison and VCF population comparison; choose a suitable panel/population for your study rather than assuming EUR is appropriate. |
| Dataset-wide inference | Inspect the logged `DECIDE` records for build, effect type, SE scale, P-value scale, frequency interpretation and strand consensus. Insufficient evidence can stop the run rather than force a decision. |
| Variant retention | Invalid records, duplicate groups, unresolved strand cases and reference-unmatched variants can be rejected before VCF export. The default does not preserve every input row. |
| Opposite-build VCF | `policies.vcf.liftover_swap: exclude` removes successfully lifted REF/ALT-swapped records. These policy exclusions are distinct from plugin lift failures; inspect both counts. |

Use the [harmonisation configuration guide](../../modules/harmonisation/configuration.md)
to review or export the actual policies. The
[processing order](../harmonisation/processing-order.md) explains each decision
and rejection boundary; do not change a threshold only to improve retention.

After a successful run, inspect:

- `results/harmonisation/STUDY/harmonisation/STUDY_harmonisation_report.html`;
- the dataset's `rejected/` records and `qc_summary/STUDY_reject_reasons.tsv`;
- the run-level summary and HTML report under
  `results/harmonisation/run_metadata/`.

Compare retained/rejected counts with the original study, read the reasons for
large losses, and inspect source-build versus target-build liftover accounting.
The two final VCFs need not contain identical numbers of variants. QC rule counts
can overlap; adding them does not give the number of distinct rejected variants.

With the packaged output layout, the GRCh37 VCF for the next steps is:

```text
results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz
```

Keep its index beside it. Only continue after the required outputs and
completion state validate. Downstream study-VCF commands require the PostGWAS
provenance declarations; an arbitrary GWAS-VCF is not interchangeable. Do not
add headers manually to bypass that boundary. See [input and output contracts](../core/input-output-contracts.md).

## 4. Review QC without changing the VCF

```console
postgwas qc \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results/qc
```

Open `results/qc/reports/STUDY_GRCh37_qc_report.html` and review
the assessed input, rule counts, missing values, and frequency/INFO provenance.
QC reports which records satisfy the configured rules but **does not filter the
input VCF**. This example therefore passes the same VCF to MAGMA below.

If your analysis requires exclusion rules, run [sumstat filtering](../modules/filtering.md)
explicitly, inspect its retention evidence, and use its validated filtered VCF
instead. Do not choose thresholds only to make a later analysis succeed.

## 5. Run MAGMA in pipeline mode

Inspect the selected plan first:

```console
postgwas pipeline --modules magma --help
```

For this target the plan is `formatter → magma`. PostGWAS generates and passes
the SNP-location and P-value/sample-size tables; you still supply both external
references:

```console
postgwas pipeline \
  --modules magma \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/GRCh37_EUR_reference \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results/magma_pipeline
```

Here `--resolve-variants-to-reference` explicitly checks formatter IDs against
the BIM identifier field and retains exact matches under the configured overlap
policy. It is not a coordinate liftover or a substitute for compatible alleles,
build, and ancestry. Review the reported overlap and exclusions before
interpreting the gene results.

For this positional example, the packaged gene window is 35 kb upstream and
10 kb downstream; the packaged MHC policy excludes both MHC SNPs and genes.
Review these choices, the tested chromosomes, and the gene-exclusion audit in
the resolved MAGMA configuration and report before interpreting association
results.

With this two-step plan and the packaged layout, inspect:

| File below `results/magma_pipeline/02_magma/` | What it tells you |
|---|---|
| `04_reports/STUDY_magma_report.html` | Run status, reference overlap, mapping and gene-analysis summaries, and exclusions |
| `03_results/positional/STUDY_magma_genes_annotated.tsv` | Tested genes, native association results and configured multiple-testing corrections |
| `02_intermediates/positional/` | Annotation, native MAGMA results and supporting execution evidence |

Read both unadjusted and corrected association results, check the number of
tested genes and the correction family, and account for excluded genes. A
significant gene association is not proof that the gene causes the trait.
There are no expected example hits because this walkthrough uses your study.
The [MAGMA guide](../modules/magma.md#outputs) gives the complete output contract.
A gene-set file is not needed for this gene-only example; competitive gene-set
analysis is a separate option with its own gene-ID compatibility requirements.

## 6. Choose the next analysis or a different execution mode

- To run exactly the same MAGMA steps yourself, use the
  [paired direct-mode recipe](../core/running-modules-independently.md#example-magma-with-and-without-pipeline-orchestration).
- To combine targets or select other engines, read
  [Pipeline Workflow](../core/pipeline-workflow.md). Dependencies depend on the
  selected method, not only the module name.
- For reusable settings, read [Configuration](../core/configuration.md). A run
  YAML can hold analysis choices and reference paths where the command supports
  it; the sample sheet still describes the raw study.
- Before archiving or rerunning, review [Output Structure](../reference/output-structure.md)
  and [Logging and Reproducibility](../core/logging-and-reproducibility.md).

If a step fails, stop at the first error and follow
[Troubleshooting](../help/troubleshooting.md). File existence alone is not a
successful or resumable analysis.
