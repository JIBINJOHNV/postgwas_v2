# GCTA fastBAT and mBAT-combo

The `gcta_gene` module runs the fastBAT and mBAT-combo association tests
implemented by GCTA. It exposes every analysis mode documented on the official
fastBAT and mBAT-combo pages:

- `fastbat_gene` tests gene regions defined by a four-column gene list.
- `fastbat_segment` tests consecutive fixed-size genomic segments.
- `fastbat_set` tests custom sets of SNPs, including pathway-defined sets.
- `mbat_combo` combines the signed mBAT and unsigned fastBAT evidence. With
  `print_component_p_values: true`, the GCTA output contains `P_mBATcombo`,
  `P_mBAT`, and `P_fastBAT` from the same variants and LD matrix. PostGWAS does
  not run a redundant second fastBAT command in this mode.

All fastBAT modes combine single-variant association p-values while accounting
for LD estimated from the configured individual-level PLINK reference.

## Input contract

All fastBAT modes and mBAT-combo consume the same GCTA `.ma` schema:

```text
SNP A1 A2 freq BETA SE P N
```

The mappings from harmonised GWAS-VCF are:

| GCTA field | PostGWAS source | Meaning |
|---|---|---|
| `SNP` | resolved variant ID | Identifier matched to the PLINK BIM file |
| `A1` | `ALT` | Effect allele |
| `A2` | `REF` | Other allele |
| `freq` | `FORMAT/AF` | Frequency of effect allele A1 |
| `BETA` | `FORMAT/ES` | Signed effect estimate |
| `SE` | `FORMAT/SE` | Standard error |
| `P` | `10^-FORMAT/LP` | Raw association p-value |
| `N` | `FORMAT/SS` | Total study sample size |

Effective sample size is not substituted for `N`. Records missing a required
value, with non-finite effects, non-positive SE or N, or with frequency outside
the open interval `(0, 1)` cannot be represented safely in GCTA `.ma` input.
The formatter reports how many such records were excluded.
Biallelic indels are retained: A1 and A2 must be nonempty and distinct, and
their complete allele strings must match the PLINK BIM pair. Strand-complement
matching, when explicitly enabled, remains limited to single-base A/C/G/T SNPs.

## Compatibility validation

The module requires one explicit `genome_build` declaration. That declaration
applies jointly to the GWAS, PLINK LD reference, and gene coordinates, which
must all use the same build. The reference population must also be explicit;
PostGWAS does not infer ancestry from filenames. Before executing GCTA it
validates:

- nonempty `.bed`, `.bim`, and `.fam` files;
- unique BIM and summary-statistic variant IDs;
- execution-mode-appropriate matching to BIM column 2;
- allele-pair compatibility for overlapping GCTA `.ma` variants;
- valid, unique four-column gene records;
- chromosome overlap between gene coordinates and the LD reference; and
- GCTA version at or above the configured minimum.

In direct mode, PostGWAS compares every `SNP` value in `--gcta-input-file`
against PLINK BIM column 2 and reports the number of unique IDs in each file,
the exact shared IDs, summary-statistic IDs absent from the BIM, and BIM IDs
absent from the summary statistics. Summary-statistic IDs absent from the BIM
are reported as a compatibility warning because GCTA will not use those variants.
PostGWAS passes the original table to GCTA without renaming IDs, removing rows,
or writing a replacement summary-statistics table. The analysis stops when no
exact IDs are shared or when a shared ID has an incompatible allele pair; a
partial exact-ID overlap alone does not stop direct mode and is not evaluated
against the pipeline reconciliation minimum.

In pipeline mode, the module uses the same harmonised GWAS-VCF supplied to the
formatter to reconcile identifier namespaces. Exact BIM IDs have priority.
When `variant_harmonisation.coordinate_fallback` is enabled, a summary variant
whose rsID is absent from the BIM can be renamed to BIM column 2 only when its
chromosome, position, and unordered REF/ALT pair identify one unambiguous BIM
record. This supports references whose BIM uses either rsIDs or coordinate-based
unique IDs without parsing or synthesizing an identifier. Direct ID matches
with conflicting coordinates or alleles stop the analysis. Ambiguous matches
remain unresolved, unresolved rows are reported and removed, and the analysis
stops if the resolved fraction is below
`variant_harmonisation.minimum_overlap_fraction`. The exact GCTA-ready table is
saved at `output_layout.harmonised_input`, and all matching counts are recorded
in the canonical log. The packaged `0.50` minimum applies only to pipeline
reconciliation and is a configurable PostGWAS fail-fast safety policy, not a
threshold defined by GCTA; the resolved fraction is always reported so a study
can adopt a stricter policy.

`fastbat_gene` and `mbat_combo` require `gene_annotation.file`.
`fastbat_segment` requires the positive configured `segment_size_kb` and does
not use an annotation file. `fastbat_set` requires exactly one set source:
`set_annotation.file`/`--fastbat-set-list` in the official block format, or
`set_annotation.gmt_file`/`--gmt` for automatic conversion:

```text
SET_1
rs123
rs456
END
```

Each set must terminate with `END`, contain unique variant IDs, and retain at
least one variant shared by the GWAS input and LD reference. Requested, matched,
unmatched, and analysis set-variant counts are logged. GCTA 1.94.1 and 1.95.3
have a compiled limit of 20,000 variants per custom set and otherwise warn and
silently ignore the set. PostGWAS checks the schema-validated
`set_annotation.maximum_set_variants` after GWAS/BIM intersection and
within-set deduplication. This GWAS-specific limit is deliberately deferred
when the standalone utility creates a reusable reference-wide set resource;
the direct or pipeline analysis applies it after intersecting that resource
with the selected GWAS. `set_annotation.oversized_set_policy` explicitly
chooses `omit` (the packaged default, with names and counts logged) or `error`.
A prepared set list bypasses GMT conversion completely. `--gmt` is rejected for every method other than
`fastbat_set`; PostGWAS does not infer or change the analysis method from the
input type.

GCTA 1.94.1 prints its version banner but returns a nonzero status for a
version-only command because no `--out` argument was supplied. PostGWAS treats
the configured banner pattern as the version-probe contract and records the
probe exit code; the actual gene-test command still requires a zero exit code
and all configured outputs.

GCTA performs its configured A1-frequency comparison for mBAT-combo. Removed
variants are retained in the official `.freq.badsnps` output when GCTA creates
that file.

## Gene-list resource preparation

Gene-list creation is deliberately outside the runtime module. The supporting
utility `tools/resource_preparation/gcta_gene_list.py` converts an explicitly
described source table into the shared headerless four-column GCTA format and
writes a provenance manifest, README, validation metrics, and SHA-256 checksums.
Keeping this one-time curation step outside `postgwas.modules.gcta_gene` prevents
an analysis run from silently changing its configured reference resource.

The source column positions, allowed chromosome labels, genome build, and output
location are required command-line inputs. The utility does not infer a genome
build, discard invalid rows, expand gene boundaries, or modify coordinates. Its
generated `resource_manifest.yaml` must travel with the gene list. Configure the
finished file through `gene_annotation.file`; the same file is used for
`fastbat_gene` and `mbat_combo`. Publication is all-or-nothing, and the utility
refuses an existing output directory so a curated resource cannot be replaced
implicitly.

## GMT pathway preparation for fastBAT

For a normal analysis, pass `--gmt pathways.gmt` with `--method fastbat_set`,
`--gene-list`, and `--gcta-reference-prefix`. PostGWAS converts the GMT file
before GCTA runs and stores the prepared set, mapping reports, manifest, and
checksums under the configured `output_layout.prepared_set_directory`. A
resumed run reuses this resource only when its source-file hashes, analyzable
variant count, genome build, gene window, conversion policies, generated-set
checksum, and every tracked audit-output checksum still match.

For GMT input, PostGWAS forms the exact analyzable universe as the validated
GWAS/BIM intersection before mapping BIM positions to genes. Pathway variants
are then deduplicated, empty and oversized policies are applied to those final
analyzable memberships, and the GCTA block file is written. The configured
`output_layout.analysis_set_list` is still independently validated from that
resource before execution. This defence-in-depth check matches GCTA's SNP
availability rules while preventing a zero-SNP or oversized set from reaching
fastBAT. Sets with no shared variants follow
`set_annotation.conversion.empty_pathway_policy`: `omit` records and excludes
them, while `error` stops before GCTA. The canonical log records input,
retained, and omitted set counts and variant-membership counts. The GMT,
gene-coordinate, BIM, and exact analysis-variant source fingerprints are
preserved in the resource manifest.

The supporting utility `tools/resource_preparation/gmt_to_gcta_fastbat_sets.py`
performs the same conversion when a reusable set resource is wanted before an
analysis. It converts standard GMT gene sets into the SNP-set block format
required by `--fastBAT-set-list`. This
is an identifier-to-variant mapping step rather than a syntax-only conversion: each GMT
gene is resolved through a build-matched four-column gene list, the configured
gene window is applied, and overlapping variants are selected from a sorted
PLINK BIM file. The output copies BIM column 2 exactly, whether it contains an
rsID or another unique variant ID; it never synthesizes or translates IDs.

Chromosome normalization and the policies for duplicate GMT genes, unmapped
genes, and pathways with no BIM variants are all explicit required arguments.
The BIM must contain unique IDs and contiguous chromosome blocks with
nondecreasing positions. Variants covered by multiple genes are written once
per pathway. The default normalized audit stores `pathway -> gene` in
`pathway_gene_mapping.parquet` and `gene -> variant` in
`gene_variant_mapping.parquet`; joining these tables reconstructs the complete
relationship without repeating every gene-variant pair for every pathway.
`conversion.audit.level: summary` omits relationship tables, while
`expanded` additionally writes `pathway_gene_variant_mapping.parquet` and must
be selected explicitly because highly overlapping pathway collections can make
that Cartesian expansion extremely large. All Parquet output is written in
bounded configured batches with configured compression. The utility also
produces pathway summaries, unmapped-gene records, a YAML manifest, and
checksums, and publishes the directory only after every validation succeeds.

For full pathway collections, conversion scans GMT one pathway at a time,
indexes eligible BIM variants into temporary chromosome caches, maps
each requested gene interval with binary searches, and stores compact integer
variant indexes with one shared BIM-ID dictionary. It then streams one pathway
at a time while writing the required GCTA block file and bounded audit batches.
An analysis-integrated conversion indexes only exact GWAS/BIM matches; the
standalone reusable-resource utility indexes the full validated BIM and defers
GWAS-specific size filtering until analysis.
This avoids retaining all GMT gene memberships and avoids repeating variant-ID
strings in every gene cache. Chromosome workers are bounded by
`execution.threads`, the resolved `execution.memory_gb`, actual
chromosome-cache sizes, and
`set_annotation.conversion.parallelism.worker_memory_multiplier`. Temporary
caches remain inside the unpublished staging directory and are removed before
atomic publication.

Before output writing, conversion calculates retained, empty, and oversized
pathway counts; exact GCTA block bytes; and a conservative uncompressed-size
estimate for configured audit outputs. It compares that estimate with the
target filesystem using `conversion.disk.minimum_free_gb` and
`conversion.disk.estimation_safety_factor`. Insufficient space stops before
large outputs are created. A later filesystem-full error reports the failed
file, staging directory, bytes written, estimate, remaining space, audit level,
and safe restart boundary; incomplete staging is never published.

With the packaged `logging.show_progress: true`, a direct `fastbat_set` run
displays seven measurable GCTA stages beneath the mandatory one-operation
module progress. Pipeline reconciliation with its GWAS-VCF displays nine. GMT input
adds eight nested preparation stages: input/GMT scanning, gene-boundary
resolution, analyzable BIM mapping, pathway-membership validation, pathway and
disk preflight, set and normalized-audit writing, provenance and checksum
writing, and atomic publication. The BIM
mapping stage has its own measured counter: its denominator is the validated
reference-variant count from the earlier PLINK stage, and its numerator advances
only after a BIM row has been fully mapped. Refreshes follow
`logging.progress_refresh_seconds`, so the counter does not rescan the BIM or
refresh the terminal for every variant. The active GCTA execution stage contains
another measured counter. For mBAT-combo,
its denominator comes from GCTA's native mapped-gene plan and its numerator
counts only complete result rows appended after successfully analysed genes.
If GCTA omits an invalid gene result, this makes the running count a
conservative lower bound rather than an invented estimate. For
fastBAT, which writes its result table only after all tests are computed, the
counter reads GCTA's native `N of total genes/sets` compute milestones; a
validated custom-set plan supplies its initial denominator. It therefore
reports measured genes, segments, or sets rather than reusing the outer stage
percentage. The counter remains below 100% until GCTA exits successfully and
PostGWAS validates the result schema. Failures preserve the last measured
count. The durable screen log records the same progress milestones without
terminal redraw characters, and
`logging.progress_refresh_seconds` controls the append-only sampling interval.

Example:

```console
python tools/resource_preparation/gmt_to_gcta_fastbat_sets.py \
  --gmt pathways.gmt \
  --gene-list NCBI37.3.gcta.gene_symbol.txt \
  --bim reference.bim \
  --output-directory prepared_pathways \
  --output-name pathways.fastbat.set \
  --resource-name "GRCh37 pathway fastBAT sets" \
  --genome-build GRCh37 \
  --gene-window-kb 50 \
  --allowed-chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
  --chromosome-label-policy exact \
  --duplicate-gene-policy deduplicate \
  --unmapped-gene-policy report \
  --empty-pathway-policy omit \
  --threads 6 \
  --memory-gb 16
```

Use a separately generated file with
`--method fastbat_set --fastbat-set-list PATH`; this direct prepared-set mode
does not run the GMT converter. It still creates the analysis-specific
intersection file described above, because prepared sets can contain BIM IDs
that are absent from a particular GWAS.
It must not be presented as an mBAT-combo pathway analysis: the official GCTA
mBAT-combo interface has no arbitrary SNP-set option corresponding to
`--fastBAT-set-list`.

## Configuration

The canonical defaults are in
`src/postgwas/config/defaults/modules/gcta_gene.yaml`. Important settings are:

- `method`: `fastbat_gene`, `fastbat_segment`, `fastbat_set`, or `mbat_combo`;
- `genome_build`, shared by all three coordinate-bearing inputs;
- `reference.prefix` and `reference.population`;
- `variant_harmonisation`, including coordinate fallback, chromosome-label,
  minimum-overlap, and strand-complement policies;
- `gene_annotation.file`;
- `set_annotation.file` or `set_annotation.gmt_file`, which are mutually exclusive;
- `set_annotation.maximum_set_variants`, the validated GCTA custom-set limit;
- `set_annotation.oversized_set_policy`, which controls sets above GCTA's
  20,000-variant hard limit;
- `set_annotation.conversion`, including chromosome, duplicate-gene,
  unmapped-gene, empty-pathway, bounded chromosome-worker, audit level,
  Parquet compression/batch size, disk preflight, and output-name policies;
- `segment_size_kb`;
- `gene_window_kb`;
- `reference_maf_min`;
- `fastbat_ld_cutoff`;
- `mbat_svd_gamma`;
- `frequency_difference_max`;
- `print_component_p_values`;
- `write_snpset`, mapped to the official mode-specific fastBAT or mBAT flag; and
- `reporting`, which controls the number of ranked associations displayed,
  p-value precision, nominal, family-wise, and FDR alpha levels, correction
  column names, and method-specific chromosome columns.

The GCTA analysis thresholds declared in the packaged YAML follow documented
GCTA defaults or examples and remain user-selectable; PostGWAS validation
thresholds are identified separately above. The exact resolved
configuration, command, GCTA version, overlap counts, allele checks, gene
counts, result columns, and generated artifacts are written to the canonical
log. CLI help displays the corresponding configured YAML default wherever an
option has one, but keeps the argparse default suppressed; therefore the YAML
remains the only source of runtime defaults.

The normalized result retains every original GCTA column and appends configured
columns for nominal significance, Bonferroni-adjusted p-values and significance,
and Benjamini-Hochberg (BH) FDR-adjusted p-values and significance. Corrections
use the complete family of tested units in that result, preserve row order, and
never filter rows. Bonferroni adjusted p-values are `min(m*p, 1)`. BH adjusted
p-values use ranked raw p-values and a reverse cumulative minimum, matching the
standard `p.adjust(method = "BH")` definition.

The completion summary reports tested units, GWAS/reference variant coverage,
the smallest p-value, and the counts meeting nominal, BH-FDR, and Bonferroni
criteria. It also lists the lowest-p-value genes, segments, or sets with raw,
BH-adjusted, and Bonferroni-adjusted p-values. For mBAT-combo, ranked entries
also show `P_mBAT` and `P_fastBAT` when GCTA was configured to emit them.
Analysis warnings identify unresolved GWAS variants and custom sets omitted
by the configured empty- or oversized-set policies. For gene, mBAT-combo, and
segment results, the reporting configuration also declares the result chromosome
column; PostGWAS reports chromosome coverage and warns when a chromosome shared
by the annotation/reference inputs has no tested units in the result. The same
findings and each ranked association are recorded in the canonical log.
Reporting-only changes do not invalidate an otherwise matching completion
manifest or rerun GCTA. If correction alpha levels or output-column names
change, a resumed run atomically regenerates the normalized table from the
checksum-validated raw GCTA result and refreshes its completion-manifest
fingerprint without re-executing GCTA.

After a successful analysis, `output_layout.completion_manifest` records the
resolved module configuration, GCTA version, SHA-256 fingerprints of the exact
summary input, annotation, PLINK reference companions, raw result, normalized
result, and validated result metrics. Resume skips GCTA only when every recorded
value still matches; a missing or changed manifest, input, or result fails with
an instruction to review and use `--overwrite`.

The configured primary outputs follow the official list: `.gene.fastbat` for
gene tests, `.seg.fastbat` for segment tests, `.fastbat` for custom sets, and
`.gene.assoc.mbat` for mBAT-combo. The same GCTA page also labels detailed
fastBAT schemas with `.fbat`; therefore every expected filename remains
version-overridable in `output_layout.primary_results` and should be protected
with a deployment-specific integration test.

## Examples

Create inputs from a harmonised GWAS-VCF:

```console
postgwas formatter \
  --vcf study.vcf.gz \
  --format gcta_gene \
  --dataset-id STUDY \
  --output-directory formatted
```

Run mBAT-combo directly:

```console
postgwas gcta_gene \
  --method mbat_combo \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --gcta-reference-prefix reference/1000G_EUR \
  --gene-list genes_grch37.txt \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Run segment-based fastBAT:

```console
postgwas gcta_gene \
  --method fastbat_segment \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --fastbat-segment-size-kb 100 \
  --gcta-reference-prefix reference/1000G_EUR \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Run custom-set fastBAT:

```console
postgwas gcta_gene \
  --method fastbat_set \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --fastbat-set-list pathway_snps.txt \
  --gcta-reference-prefix reference/1000G_EUR \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Run pathway-set fastBAT directly from GMT:

```console
postgwas gcta_gene \
  --method fastbat_set \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --gmt pathways.gmt \
  --gene-list genes_grch37.txt \
  --gcta-reference-prefix reference/1000G_EUR \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Run through the PostGWAS pipeline:

```console
postgwas pipeline \
  --modules gcta_gene \
  --method mbat_combo \
  --vcf study.vcf.gz \
  --gcta-reference-prefix reference/1000G_EUR \
  --gene-list genes_grch37.txt \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Pipeline resume is enabled by default. When the formatter step completed in an
earlier attempt, PostGWAS verifies its completion manifest, input VCF checksum,
resolved formatter configuration, output schemas, and output checksums before
supplying the saved GCTA input to this module. Existing GCTA and GMT-conversion
results are likewise reused only through their validated resume contracts. Pass
`--no-resume` to disable continuation for a run.

## References

- [Official GCTA fastBAT documentation](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT)
- [Official GCTA mBAT-combo documentation](https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo)
- [Official GCTA mBAT-combo source](https://github.com/JianYang-Lab/GCTA/blob/main/main/mbat.cpp),
  which logs the mapped-gene total and flushes one output row after each tested
  gene;
- [Official GCTA fastBAT source (`sbat.cpp`)](https://github.com/JianYang-Lab/GCTA/blob/main/main/sbat.cpp),
  which resolves summary and set entries through the reference SNP-name map;
- Bakshi et al. 2016, [doi:10.1038/srep32894](https://doi.org/10.1038/srep32894)
- Li et al. 2023, *American Journal of Human Genetics* 110:30–43,
  [doi:10.1016/j.ajhg.2022.12.006](https://doi.org/10.1016/j.ajhg.2022.12.006)
- Benjamini and Hochberg 1995, *Journal of the Royal Statistical Society B*
  57:289–300,
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [R `p.adjust` reference implementation and definitions](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
