# MAGMA gene and gene-set analysis

PostGWAS runs single-trait MAGMA gene association and, when a gene-set file is
provided, competitive gene-set analysis. The canonical source of runtime
values is `src/postgwas/config/defaults/modules/magma.yaml`; an exported module
or pipeline YAML can override those values, and explicit CLI arguments override
the corresponding YAML keys.

## Inputs and scientific decisions

The formatter creates two headed audit tables. The p-value table contains the
configured variant ID, raw p-value, and per-variant total sample size. The
SNP-location table retains chromosome, position, `REF`, and `ALT` so reference
compatibility can be checked without reopening the GWAS-VCF. Immediately before
`--annotate`, the MAGMA runner emits the required headerless three-column
`SNP CHR BP` file.
The whitespace character used for both prepared MAGMA tables comes from
`input.output_table_delimiter`; schema validation limits it to a tab or space,
which are accepted by MAGMA's whitespace-delimited input contract.

Variant-ID selection is a general formatter feature, not MAGMA-specific code:

```yaml
variant_identifiers:
  default_type: rsid
  target_types: {}
  rsid_pattern: '(?i)^rs[0-9]+$'
  rsid_extraction_pattern: '(?i)(?:^|;)(rs[0-9]+)(?:;|$)'
  unique_id_template: '{chromosome}_{position}_{reference_allele}_{alternate_allele}'
```

`rsid` extracts one configured rsID from the VCF `ID` field. `unique` constructs
the configured chromosome-position-REF-ALT identifier from validated VCF
fields. `target_types` may choose a different convention for each formatter
target, so selecting MAGMA-compatible IDs does not alter LDSC, MiXeR, or another
export. Standalone formatter runs can override the common choice with
`--variant-id-type {rsid,unique}`.

In a MAGMA pipeline, PostGWAS scans every row of BIM field 2 before formatting.
A homogeneous rsID reference selects `rsid` for the MAGMA target; a homogeneous
configured coordinate-and-allele convention selects `unique`. Mixed or custom
identifier systems fail with their observed counts instead of guessing. PLINK
BIM fields 5 and 6 are A1/A2 rather than guaranteed REF/ALT, so classification
recognizes either allele order. The formatter itself still constructs the
configured VCF REF/ALT order, and optional exact intersection later confirms
that the resulting IDs, coordinates, and unordered allele pairs agree.

Reference-ID resolution is opt-in. Its canonical default is:

```yaml
snp_harmonisation:
  resolve_variants_to_reference: false
```

With this default, PostGWAS validates and deduplicates the formatter tables but
does not pre-filter them against the BIM. MAGMA receives those prepared files
directly and performs its normal reference-based analysis. Therefore
`minimum_overlap_fraction` is not applied in this mode. This is appropriate
when the formatter has already produced the intended reference-ID convention.

Exact pre-intersection can be enabled in YAML or from the CLI:

```console
postgwas magma ... --resolve-variants-to-reference
```

When enabled, PostGWAS retains only formatter variants whose ID occurs in BIM
field 2 and whose normalized chromosome, integer position, and unordered allele
pair agree with that BIM record. An ID-only match with a different coordinate
or allele pair is excluded and reported separately; it is not silently treated
as a strand or coordinate match. The compatible unique-variant fraction must
meet `minimum_overlap_fraction`; otherwise analysis stops with the observed and
required fractions and the mismatch counts. This operation is deliberately an
intersection only: it does not construct IDs, read synonym files, or rescue an
incompatible record.

In either mode, duplicate formatter IDs are consolidated by the configured
`lowest_p` policy and counted. MAGMA is then run with `duplicate=error`, so a
remaining duplicate cannot be silently discarded by the external tool.

The configured overlap fraction is a PostGWAS safety threshold, not a MAGMA
statistical parameter. It prevents an analysis using a largely incompatible LD
reference. Before gene association, PostGWAS checks gene-set identifiers against
the complete gene universe of the selected mapping reference. It does not use the
smaller, study-specific annotation produced after intersecting variants with that
study, because genomic coverage is not evidence of an identifier-system mismatch.
After gene association, PostGWAS separately records how many tested genes occur in
the supplied gene sets and how many sets contain at least one tested gene. This
study-specific coverage is descriptive and is not subjected to the identifier
compatibility threshold; zero represented sets remains a fatal input error.

This decision follows the official
[PLINK BIM format specification](https://www.cog-genomics.org/plink/1.9/formats#bim),
which defines field 2 as the variant identifier and fields 5/6 as A1/A2; the
[VCF 4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf),
which defines semicolon-separated identifiers in the `ID` field; and the
official [MAGMA documentation and resources](https://cncr.nl/research/magma/).
This separation also follows MAGMA's competitive gene-set model, which compares
genes in a set with the other genes available in the gene-analysis results, as
described in the
[primary MAGMA publication](https://doi.org/10.1371/journal.pcbi.1004219).

Gene sets may use standard GMT records (`name`, `description`, then tab-separated
gene IDs), native MAGMA set-annotation records (`name`, then space-separated
gene IDs), or a configured two-column set/gene membership table.
`gene_sets.input_format` can require a format or use the logged `auto` detection
for GMT/native inputs. PostGWAS writes one validated native MAGMA set file before
calling the external tool. Native MAGMA and membership files do not carry
descriptions, so
`gene_set_description` is null for those records. The output column names for
the description and supplied genes come from `result_schema`; the original
MAGMA `VARIABLE` column is retained, and `FULL_NAME` is retained when present
or added as a copy of `VARIABLE` when MAGMA omits it.

`N_COL` is total sample size, including for case-control GWAS. The official
MAGMA manual recommends `ncol` when per-variant sample size is available,
particularly for meta-analysis. Effective sample size is not substituted.

Supported gene-model configuration is limited to models documented for SNP
p-value input: `snp-wise`, `snp-wise=mean`, `snp-wise=top`,
`snp-wise=top,P`, `multi=snp-wise`, `snp-wise=multi`, and `multi`. Models that
require raw genotype/phenotype analysis are rejected.

Genome build and population are currently recorded in configuration and logs
and each selected mapping definition must declare the same values. Because BIM
and MAGMA annotation files do not provide authoritative build/population
metadata, users must still ensure the files themselves match those declarations.

## Functional mapping methods

`mapping.selected` accepts any configured mapping-definition names, and
`mapping.primary` chooses the one passed to downstream modules. Each selected
mapping is run independently and
written below `mappings/<mapping-name>/`. The combined mapping table is a
provenance-rich result catalogue; p-values from different annotations are not
pooled or treated as exchangeable.

- `positional` creates the ordinary MAGMA SNP-to-gene annotation from the
  configured gene locations and windows.
- `emagma` consumes a tissue-specific eQTL SNP-to-gene annotation. Optional
  tissue network memberships can drive competitive gene-set testing.
- `h_magma` consumes a tissue/cell-specific chromatin-interaction annotation.
- `n_magma` creates the positional component for the current study, then unions
  it with all configured Hi-C, eQTL, and TOM co-expression components. Repeated
  SNP assignments within a gene are removed. The configured nMAGMA
  `gene_location_file` supplies the authoritative protein-coding gene coordinates;
  component coordinate differences, invalid component coordinate placeholders,
  and genes absent from that reference are counted in the log. The packaged
  nMAGMA definitions use a zero-kilobase
  positional window, matching upstream nMAGMA rather than inheriting the
  conventional MAGMA window.
- `chrom_magma` annotates and tests regulatory elements, links element results
  to genes, and selects the lowest-p element per gene with a deterministic
  element-ID tie break. This published gene mapping is reported separately and
  is not relabelled as an ordinary MAGMA gene-body result. Selecting the minimum
  of multiple element p-values does not produce a calibrated gene p-value, so
  PostGWAS does not apply ordinary gene-level Bonferroni/FDR correction or report
  gene significance for this ranking. Competitive MAGMA gene-set testing is
  therefore rejected for a chromMAGMA definition.

Choose a positional, eMAGMA, H-MAGMA, or nMAGMA definition as `mapping.primary`
when a later pipeline module consumes ordinary MAGMA gene units. chromMAGMA may
run in the same invocation and remains available in its mapping-specific output.

External eMAGMA/H-MAGMA and merged nMAGMA annotations are structurally validated
and must meet `annotation_validation.minimum_bim_variant_overlap_fraction` for
exact SNP IDs in BIM field 2. This is independent of the optional formatter
intersection controlled by `snp_harmonisation.resolve_variants_to_reference`.
The pinned bundle records the identifiers actually present upstream: Entrez IDs
for eMAGMA, Ensembl IDs for H-MAGMA, gene symbols for nMAGMA, and mixed target
identifiers for chromMAGMA. PostGWAS preserves and reports these types; it does
not silently relabel functional-mapping targets.

## Execution and outputs

MAGMA v1.10 or newer is required and is identified using the official
`--version` interface. Commands, version, resolved configuration,
input compatibility counts, duplicate resolution, batching decisions, output
paths, and failure details are written to the canonical log. Batch count is
bounded by configured threads, total memory, memory per MAGMA process, and
minimum genes per batch.

The terminal shows one live overall progress bar driven by the same stage
contexts used by the canonical log. It names the active stage, retains a timed
completion line and concise scientific outcome for each finished stage, and
identifies the failed stage if the run stops. Variant outcomes distinguish
optional exact-BIM exclusions from duplicate rows consolidated by the configured
lowest-p-value rule. Gene and gene-set outcomes report the tested count,
nominal significance, and each configured global correction when the mapping
produces calibrated MAGMA gene p-values. chromMAGMA instead reports the number
of mapped genes and its ranking-statistic interpretation. The reporting
cutoff comes from
`multiple_testing.reporting_significance_threshold`; it does not filter or
otherwise change a result. The correction order follows the configured method
lists, and their terminal names come from
`multiple_testing.reporting_method_labels`. Outcome values are displayed as
indented key-value fields with one aligned colon column across every MAGMA
stage; wrapped values continue under the value column rather than returning to
the left margin. The shared global `logging.terminal_label_width` setting
controls that label column for progress and completion summaries. Set
`logging.show_progress: false` in a full run configuration to suppress this
terminal-only display; the structured stage outcomes are always retained in
the canonical log.

All analysis files are first written below the configured `.partial` staging
directory. They are published only after every required output has been
validated. A failure leaves the isolated partial directory for diagnosis and
finalizes the canonical log; it does not expose a partial gene result as a
completed result.

The human-readable MAGMA gene table is retained unchanged as the raw
`.genes.out` result. PostGWAS also writes the configured
`corrected_genes` TSV, preserving every gene row and original column and adding
Bonferroni and Benjamini–Hochberg adjusted p-values across the complete family
of valid gene p-values. Missing, non-finite, or out-of-range gene p-values fail
validation so the multiple-testing denominator cannot change silently. The
methods and appended column names come from `multiple_testing.gene_methods` and
`result_schema.global_correction_column_pattern`. Here, FDR means the standard
Benjamini–Hochberg procedure (`fdr_bh`), not a local-FDR estimate. The official
MAGMA manual defines `P` as the gene p-value in the default `.genes.out` schema;
`result_schema.gene_p_value_column` remains configurable for another supported
MAGMA result schema.

The mapping catalogue adds configured columns identifying the mapping name,
method, biological context, gene-ID system, annotation source/version, statistic
type, and statistical interpretation. This prevents chromMAGMA's minimum linked
element p-value from being confused with a calibrated MAGMA gene p-value.

Lines beginning with `#` in MAGMA's human-readable result files are tool
metadata rather than result rows. PostGWAS treats that prefix as part of the
MAGMA output contract, just like `.genes.out` and `.gsa.out`, rather than as a
runtime analysis parameter.

Gene-set reports contain configured global and named-family Bonferroni,
Šidák, Holm, and Benjamini–Hochberg adjustments. Each tested family and its
denominator are logged. Named family regular expressions and correction methods
come only from YAML.

## Configuration and commands

```sh
postgwas config export \
  --module magma \
  --style full \
  --output magma.yaml

postgwas magma \
  --run-config magma.yaml
```

Prepare and use the pinned functional resources:

```sh
postgwas resources prepare magma

postgwas magma \
  --run-config ~/Documents/software_resources/resourses/postgwas/magma/functional_mapping/configs/magma_functional_mapping.yaml \
  --magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac \
  --primary-magma-mapping positional \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --dataset-id STUDY \
  --output-directory results
```

The preparer installs 88 schema-validated definitions: 1 positional, 47 eMAGMA,
35 H-MAGMA, 4 nMAGMA, and 1 chromMAGMA. It also installs the official GRCh37
MAGMA gene locations, dbSNP151 synonyms, and `g1000_eur` BED/BIM/FAM LD panel.
Every archive and installed file is checksum recorded in
`resource_manifest.yaml`.
The generated YAML stores these values under `modules.magma`, so the same file
is valid for both `postgwas magma` and `postgwas pipeline`.

The supplied eMAGMA annotations use Entrez IDs, whereas the accompanying GTEx
network memberships use Ensembl IDs. The preparer uses the pinned GENCODE 26
Ensembl-to-HGNC bridge and MAGMA NCBI37.3 HGNC-to-Entrez reference, retaining
only unambiguous links. Per-network mapping and annotation-overlap counts are
written to `emagma/networks_entrez/identifier_harmonisation_summary.tsv`; the
exact crosswalk is retained under `base/gene_identifiers/`. Original network
files remain unchanged under `emagma/networks/`.

Pipeline mode obtains the two summary-statistic inputs from formatter:

```sh
postgwas config export \
  --pipeline magma \
  --style full \
  --output magma_pipeline.yaml

postgwas pipeline \
  --modules magma \
  --vcf study.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config magma_pipeline.yaml
```

## Installation and licence

The Dockerfile downloads the official static Linux MAGMA v1.10 archive during
a local build, verifies its pinned SHA-256 digest, installs only the executable,
and checks its reported version. CNCR states that MAGMA versions after v1.0 use
standard copyright and that their binaries and source may not be distributed
or modified. Therefore, do not publish or redistribute an image containing the
installed v1.10 binary without permission from the MAGMA authors.

## Primary sources

- [Official MAGMA page, v1.10 downloads, manual, reference data, and licence notice](https://cncr.nl/research/magma/)
- [Official PLINK BIM format specification](https://www.cog-genomics.org/plink/1.9/formats#bim)
- [VCF 4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [de Leeuw et al. 2015, MAGMA: generalized gene-set analysis of GWAS data](https://doi.org/10.1371/journal.pcbi.1004219)
- [Gerring et al. 2021, E-MAGMA](https://doi.org/10.1093/bioinformatics/btab115)
- [Sey et al. 2020, H-MAGMA](https://doi.org/10.1038/s41593-020-0603-0)
- [Yang et al. 2021, nMAGMA](https://doi.org/10.1093/bib/bbaa298)
- [Nameki et al. 2022, chromMAGMA](https://doi.org/10.26508/lsa.202201446)
- [Benjamini and Hochberg 1995, the original false-discovery-rate procedure](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [Statsmodels reference implementation for Bonferroni and Benjamini–Hochberg adjusted p-values](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html)
