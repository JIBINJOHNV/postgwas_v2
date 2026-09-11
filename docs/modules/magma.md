# MAGMA gene and gene-set analysis

PostGWAS runs single-trait MAGMA gene association and, when a gene-set file is
provided, competitive gene-set analysis. The canonical source of runtime
values is `src/postgwas/config/defaults/modules/magma.yaml`; an exported module
or pipeline YAML can override those values, and explicit CLI arguments override
the corresponding YAML keys.

## Inputs and scientific decisions

In pipeline mode, the formatter creates the final two MAGMA tables in one
in-memory export. The p-value table contains the configured variant ID, raw
p-value, and per-variant total sample size. The SNP-location table is the
headerless three-column `SNP CHR BP` file required by MAGMA `--annotate`.
The same export measures BIM identifier overlap, optionally intersects with the
BIM, applies configured chromosome and MHC policies, and writes the paired
tables once. The MAGMA runner consumes those exact files without reopening or
rewriting them. A standalone formatter export remains a headed audit artifact
with `REF` and `ALT`; direct MAGMA validates independently supplied tables before
creating its run-owned final inputs.
The whitespace character used for both prepared MAGMA tables comes from
`input.output_table_delimiter`; schema validation limits it to a tab or space,
which are accepted by MAGMA's whitespace-delimited input contract.
Input and native-result readers use the configured
`input.table_delimiter_pattern` and `result_schema.table_delimiter_pattern`.
PostGWAS retains pandas' fast C parser for `\s+` and literal one-character
separators. A different multi-character regular expression automatically uses
the Python parser so every schema-valid configured pattern is honoured. The
selected parser, delimiter pattern, file role, and path are written to the
canonical log. This follows the
[pandas separator contract](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html),
which reserves regular-expression separators other than `\s+` for its Python
engine.

Variant-ID selection is a general formatter feature, not MAGMA-specific code:

```yaml
variant_identifiers:
  default_type: rsid
  target_types: {}
  default_duplicate_policy: exclude_all
  target_duplicate_policies: {}
  duplicate_rank_tolerance: 1.0e-12
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
They can override the duplicate fallback for every selected target with
`--duplicate-id-policy {exclude_all,error,lowest_p,most_significant,highest_maf,highest_info}`.
MAGMA uses the same shared policy as every formatter target; P-value ranking is
available only when explicitly requested and is not the packaged default.

In a MAGMA pipeline, PostGWAS scans every row of the configured BIM variant-ID
field (field 2 in the standard PLINK BIM layout) before formatting.
A homogeneous rsID reference selects `rsid` for the MAGMA target; a homogeneous
configured coordinate-and-allele convention selects `unique`. Mixed or custom
identifier systems fail with their observed counts instead of guessing. PLINK
BIM fields 5 and 6 are A1/A2 rather than guaranteed REF/ALT, so classification
recognizes either allele order. The formatter itself still constructs the
configured VCF REF/ALT order. The MAGMA compatibility gate compares the final
GWAS identifier only with the configured BIM variant-ID field; it does not reinterpret an identifier
match using BIM coordinates or alleles.

Reference-ID resolution is opt-in. Its canonical default is:

```yaml
snp_harmonisation:
  resolve_variants_to_reference: false
  duplicate_policy: lowest_p
```

With this default, PostGWAS validates and deduplicates the formatter tables but
does not pre-filter them against the BIM. MAGMA receives those prepared files
directly and performs its normal reference-based analysis. PostGWAS still
measures identifier overlap with the configured BIM variant-ID field and applies
`minimum_overlap_fraction`; the option controls filtering, not compatibility
assessment.

Exact pre-intersection can be enabled in YAML or from the CLI:

```console
postgwas magma ... --resolve-variants-to-reference
```

When enabled, PostGWAS retains only formatter variants whose ID occurs in BIM
field 2. The unique identifier-overlap fraction must meet
`minimum_overlap_fraction` in both filtering modes; otherwise analysis stops
with the observed and required fractions. This operation is deliberately an
identifier intersection only: it does not construct IDs, compare coordinates
or alleles, read synonym files, or rescue an incompatible identifier.

In pipeline mode, the formatter applies the resolved MAGMA duplicate policy
before either final table is written. The paired location and p-value tables
always contain the same resolved unique IDs in the same order. For independently
prepared direct-module inputs, the MAGMA module applies the same policy after
optional reference filtering. The accepted values are:

- `lowest_p` (default): retain one row per duplicated SNP ID using the lowest
  valid p-value; original input order resolves an exact p-value tie.
- `remove`: exclude every row belonging to each duplicated SNP-ID group.
- `err`: stop before writing the prepared MAGMA inputs, reporting the duplicated
  group and row counts plus example IDs.

Coordinate or allele conflicts in the paired location table always stop the run,
regardless of this policy. Formatter-produced inputs are normally already unique,
so the policy primarily protects independently prepared direct-module inputs.
MAGMA is always run with `duplicate=error` as a final safeguard. This prevents
MAGMA's default `duplicate=drop` behavior from automatically removing a duplicate
that remains after PostGWAS preparation. Select the policy from the CLI with
`--duplicate-policy {err,lowest_p,remove}` or set it in YAML.

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
The terminal and canonical screen log do not print a separate preflight block.
All validated inputs, quantitative identifier decisions, transformations,
analyses and outputs appear once in the canonical ordered execution summary.
Long input paths and complete schemas remain in the canonical log and HTML
report. There is no pathway-preparation stage: the validated source pathway file
is passed directly to MAGMA without modification. There is no separate
on-screen comparison of the formatter's SNP-location and p-value tables: both
are produced from the same validated retained frame, while their individual
schemas and row counts remain validated and logged.

After successful analysis, the same ordered 12-stage execution summary is
printed on screen and written to
`04_reports/<dataset>_magma_pipeline_summary.csv` and
`04_reports/<dataset>_magma_report.html`. The HTML is a scientific results
report rather than a copy of the execution table. It presents variant/reference
coverage, exclusions, annotated-versus-tested gene accounting, pathway-gene
coverage, analysis parameters, multiple-testing results, and every gene and
pathway result across searchable, sortable pages. Complete corrected and
annotated TSV files remain the authoritative result tables and are linked from
the report; large pathway membership lists are not duplicated into HTML. The
workflow audit remains available in a collapsible section. The page size and
displayed columns are controlled by `html_report` in the canonical MAGMA YAML.
The final terminal summary also lists the strongest gene and competitive
pathway associations, ranked by their unadjusted MAGMA p-values. Gene rows show
the configured gene-level corrections; pathway rows show the corrections shared
with the gene results plus the configured primary pathway correction.
`screen_summary.top_gene_rows` and `screen_summary.top_pathway_rows` control
the two list lengths; both default to 10. The ranked rows use the configured
method labels without repeating “adjusted P”, and
`screen_summary.p_value_significant_digits` controls their compact terminal
precision. Full-precision values remain in the result tables and HTML report.
These rankings are descriptive and do not replace the configured
multiple-testing significance decisions.
The CSV column names, delimiter, null marker, and both report paths come from
that YAML. Its schema version is included in the checkpointed module
configuration, so an obsolete saved screen/CSV/HTML contract cannot be reused
after that version changes. The current 12-stage report contract is schema
version 5.
VCF extraction, configured column transformations, BIM overlap, analysis-scope
policies, and final MAGMA file creation are reported together as one
variant-input preparation stage. The formatter still reads the VCF only once.

The fixed execution order is:

1. validate the input GWAS-VCF;
2. validate PLINK LD-reference files and determine the BIM identifier format;
3. validate the gene-location file;
4. validate the original pathway GMT file;
5. select compatible gene-location identifiers for pathway analysis;
6. prepare the variant-level MAGMA inputs, including VCF field transformations,
   BIM overlap, optional intersection, scope policies, and final file creation;
7. create the MAGMA SNP-to-gene annotation;
8. calculate MAGMA gene-association statistics;
9. annotate gene results and adjust gene p-values;
10. run competitive pathway analysis with the original pathway file;
11. annotate pathway results and adjust pathway p-values; and
12. validate and publish all outputs.

`gene_sets.identifier_mismatch_action` controls only optional competitive
gene-set analysis. Its default, `skip`, records a warning when the complete
mapping reference and gene-set identifiers fall below the configured overlap,
omits pathway testing, and continues through gene association and gene p-value
correction. `error` instead stops the complete MAGMA run. Malformed files and
invalid gene-association inputs remain fatal under either policy. The same choice
is available as `--gene-set-identifier-mismatch {skip,error}`.

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

## Chromosome and MHC analysis scope

MAGMA applies a local, explicit scope after optional BIM intersection and before
SNP-to-unit annotation. The packaged policy is:

```yaml
mhc:
  policy: exclude_both
  region_override: null
chromosomes:
  exclude: [Y, MT]
```

Chromosome X is retained. `exclude_snps` removes variants in the configured MHC
interval before annotation and gene testing; `exclude_genes` removes annotation
units whose tested interval overlaps the MHC before gene analysis and therefore
also removes them from the competitive-test background; `exclude_both` applies
both; `include` applies neither MHC exclusion. Every chromosome listed in
`chromosomes.exclude` is removed at both levels.

The default one-based inclusive intervals are stored under
`resources.genomes.<build>.regions.mhc`, using the Genome Reference Consortium
GRCh37/GRCh38 definitions. A run can override the complete region with
`--mhc-chrom`, `--mhc-start`, and `--mhc-end`; supplying only part of the
override fails configuration validation. Source VCF, BIM, gene-location,
precomputed `.genes.annot`, and GMT files remain unchanged. MAGMA receives a
run-owned scoped annotation. The pathway file is passed directly to MAGMA;
members absent from the gene-analysis results are discarded by MAGMA itself.

Gene sets may use standard GMT records (`name`, `description`, then tab-separated
gene IDs), native MAGMA set-annotation records (`name`, then space-separated
gene IDs), or a configured two-column set/gene membership table.
For positional and nMAGMA mappings, PostGWAS validates the configured gene-location
record as primary gene ID, chromosome, integer start, integer end, and `+`/`-`
strand, followed by an optional alternate gene ID. It compares gene-set members
with the primary ID first and does not relabel a compatible file. For positional
MAGMA, when column one is incompatible but column six reaches the configured
threshold, PostGWAS leaves every pathway identifier unchanged and creates a
run-owned location file with column six as the primary MAGMA ID and the original
column one retained as its alternate ID. Repeated column-six identifiers follow
`gene_sets.alternate_id_duplicate_policy`; the default `longest_interval` keeps
the longest interval with a deterministic primary-ID tie-break, while `error`
stops preparation. `input.alternate_gene_id_type` records the identifier system
used by the derived reference. Missing alternate IDs are omitted and all counts
are logged.
PostGWAS never substitutes, expands, or otherwise rewrites identifiers in the
pathway membership file. Non-positional mappings must therefore match the gene
identifiers produced by their supplied annotation; an incompatible pathway file
is skipped or rejected according to `gene_sets.identifier_mismatch_action`.
Analysis stops or skips pathway testing according to the configured mismatch
action when both column-one and column-six overlap remain
below the configured threshold. The overlap denominator is the smaller of the
reference and pathway gene-set universes, so a comprehensive GMT is not
penalised for containing genes absent from a protein-coding location reference.
This keeps pathway memberships unchanged while ensuring that positional MAGMA
gene results use the same identifiers. MAGMA itself
documents four required location columns and an optional fifth strand column in
the [MAGMA manual](https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf).
`gene_sets.input_format` can require a format or use the logged `auto` detection
for GMT/native inputs. PostGWAS passes the validated source file directly to
MAGMA. Headerless column-based membership files use MAGMA's native `col=`
modifier; headered membership files are rejected because MAGMA does not accept
a header for this input and PostGWAS does not rewrite pathway data. Native MAGMA
and membership files do not carry descriptions, so
`gene_set_description` is null for those records. The output column names for
the description and supplied genes come from `result_schema`; the original
MAGMA `VARIABLE` column is retained, and `FULL_NAME` is retained when present
or added as a copy of `VARIABLE` when MAGMA omits it.
This follows the MAGMA `--set-annot` contract: identifiers absent from the gene
analysis are discarded by MAGMA, and non-gene GMT fields such as description
links may remain when they cannot match the gene identifier nomenclature.

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
mapping is run independently and written below mapping-specific subdirectories
in `02_intermediates/` and `03_results/`. The combined mapping table is a
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
minimum genes per batch. Every batch uses MAGMA's native
`--batch <number> <total>` interface with identical scientific arguments and the
configured seed; MAGMA's `--merge` creates the final `.genes.raw` and
`.genes.out`. Competitive gene-set analysis consumes that merged `.genes.raw`,
which retains the gene–gene correlation information needed by MAGMA.

`--magma-memory-per-worker-gb` overrides the canonical
`modules.magma.batching.memory_per_process_gb` setting. PostGWAS chooses no more
workers than both `--threads` and `floor(--memory-gb /
--magma-memory-per-worker-gb)` permit, subject to the configured minimum genes
per batch. This value is a scheduler budget, not an operating-system memory
limit; setting it below MAGMA's actual peak use can cause out-of-memory failure.

The exclusion audit consists of a shared variant table under
`02_intermediates/exclusions/`, a scoped annotation and excluded-unit table per
mapping under `02_intermediates/<mapping>/`, and a per-mapping exclusion summary
under `03_results/<mapping>/`. The mapping comparison table repeats the MHC
policy, interval, excluded chromosomes, and excluded-unit count so the analysis
scope remains visible beside gene results.

The terminal shows one live overall progress bar driven by the same stage
contexts used by the canonical log. It names the active stage, retains a timed
completion line and concise scientific outcome for each finished stage, and
identifies the failed stage if the run stops. Variant outcomes distinguish
optional exact-BIM exclusions from the configured `err`, `lowest_p`, or `remove`
duplicate decision. Gene and gene-set outcomes report the tested count,
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
validated. The SNP-location and p-value inputs, LD-reference companions,
executable/version, and configured gene-set files are validated before that
staging directory is created. A preflight failure therefore finalizes the
canonical log without leaving a false partial-run marker. If a later run finds
a stale staging tree containing only empty directories, it records and removes
that tree, then starts normally without `--overwrite`. A failure after analysis
output has started leaves the non-empty isolated partial directory for
diagnosis and requires `--overwrite`; it does not expose a partial gene result
as completed.

The human-readable MAGMA gene table is retained unchanged as the raw
`.genes.out` result beside its paired `.genes.raw` file under the selected
mapping's `native_outputs/` directory. Keeping this pair together preserves
MAGMA's prefix-based downstream interface. PostGWAS also writes the configured
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

### Output structure

The numbered directories are defined only by `output_layout` in the canonical
MAGMA YAML:

```text
02_magma/
├── 00_run_metadata/
├── 01_inputs/
├── 02_intermediates/
│   └── <mapping_name>/
│       ├── gene_reference/
│       ├── annotations/
│       ├── batches/
│       └── native_outputs/
├── 03_results/
│   └── <mapping_name>/
├── 04_comparisons/
├── 04_reports/
├── 05_logs/
└── .partial/
```

`01_inputs` contains the validated paired MAGMA tables. Tool-native annotations,
batch files, `.genes.raw`/`.genes.out`, and `.gsa.out` remain traceable under
`02_intermediates`. Corrected and annotated tables intended for interpretation
are under `03_results`; the cross-mapping catalogue is under `04_comparisons`.
For gene-location-backed mappings, the primary interpreted gene table is
`<dataset>_magma_genes_annotated.tsv`. It contains the complete native MAGMA
statistics, the configured reference chromosome/start/end/strand and alternate
gene identifier, followed by the configured adjusted p-values. The native
MAGMA `START` and `STOP` remain the tested annotation-window coordinates; the
separate `GENE_REFERENCE_START` and `GENE_REFERENCE_END` columns are the source
gene interval.
The resolved configuration and checksum-validated completion manifest are under
`00_run_metadata`; the canonical log is under `05_logs`. Resume requires that
manifest's configuration, input, software, and output fingerprints to match.

Lines beginning with `#` in MAGMA's human-readable result files are tool
metadata rather than result rows. PostGWAS treats that prefix as part of the
MAGMA output contract, just like `.genes.out` and `.gsa.out`, rather than as a
runtime analysis parameter.

Gene-set reports contain the configured global corrections—Bonferroni, Šidák,
Holm, and Benjamini–Hochberg by default—and configured named-family corrections.
`multiple_testing.primary_method` selects exactly one global correction for the
reported significance decision; the packaged value is `fdr_bh`. Every row
records that method, its adjusted p-value, and the resulting Boolean decision in
the three configured `result_schema.primary_*` columns. Terminal stage findings
declare significance only from this primary correction. Other global and named-
family values remain in the tables as explicitly secondary sensitivity results;
they must not be selected after inspecting which method is significant.
Ordinary gene reports use the separate `multiple_testing.gene_methods` list,
which defaults to Bonferroni and Benjamini–Hochberg. Each tested family and its
denominator are logged. Correction methods and named-family regular expressions
come only from YAML. Šidák uses the stable
`-expm1(number_of_tests * log1p(-p))` calculation.

The primary correction is calculated across all tested gene sets within one
mapping analysis. The default BH-FDR interpretation retains the assumptions of
that procedure; overlapping gene sets can be dependent, so the additional
global corrections remain useful pre-specified sensitivity results. The primary
correction does not combine or correct across several selected mapping
definitions.
`mapping.primary` identifies the confirmatory top-level mapping; additional
mapping analyses must be interpreted as separately pre-specified or exploratory
unless the study defines an external cross-mapping multiplicity procedure.

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
- [Genome Reference Consortium GRCh37 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [Gerring et al. 2021, E-MAGMA](https://doi.org/10.1093/bioinformatics/btab115)
- [Sey et al. 2020, H-MAGMA](https://doi.org/10.1038/s41593-020-0603-0)
- [Yang et al. 2021, nMAGMA](https://doi.org/10.1093/bib/bbaa298)
- [Nameki et al. 2022, chromMAGMA](https://doi.org/10.26508/lsa.202201446)
- [Benjamini and Hochberg 1995, the original false-discovery-rate procedure](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [Benjamini and Bogomolov 2014, selective inference across multiple hypothesis families](https://doi.org/10.1111/rssb.12028)
- [Statsmodels reference implementation for Bonferroni and Benjamini–Hochberg adjusted p-values](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html)
