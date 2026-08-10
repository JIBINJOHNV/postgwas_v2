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

## Scientific input contract

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
- exact-ID or unambiguous coordinate-and-allele resolution to BIM column 2;
- reference overlap at or above the configured minimum;
- allele-pair compatibility for overlapping GCTA `.ma` variants;
- valid, unique four-column gene records;
- chromosome overlap between gene coordinates and the LD reference; and
- GCTA version at or above the configured minimum.

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
in the canonical log. Standalone runs should also pass `--vcf` when the
formatter input and BIM may use different identifier namespaces.
The packaged `0.50` minimum is a configurable PostGWAS fail-fast safety policy,
not a threshold defined by GCTA; the resolved fraction is always reported so a
study can adopt a stricter policy.

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
silently ignore the set. PostGWAS detects that condition before execution;
`set_annotation.oversized_set_policy` explicitly chooses `omit` (the packaged
default, with names and counts logged) or `error`. A prepared set list bypasses GMT
conversion completely. `--gmt` is rejected for every method other than
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
an analysis run from silently changing its scientific reference resource.

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
resumed run reuses this resource only when its source-file hashes, genome
build, gene window, conversion policies, and generated-set checksum still
match.

Before analysis, PostGWAS writes the configured
`output_layout.analysis_set_list` by retaining only variant IDs present in
both the harmonised GWAS input and PLINK BIM reference. This matches GCTA's
own two-stage SNP availability rules while preventing a zero-SNP set from
reaching fastBAT. Sets with no shared variants follow
`set_annotation.conversion.empty_pathway_policy`: `omit` records and excludes
them, while `error` stops before GCTA. The canonical log records input,
retained, and omitted set counts and variant-membership counts. The original
GMT-derived resource remains unchanged for provenance.

The supporting utility `tools/resource_preparation/gmt_to_gcta_fastbat_sets.py`
performs the same conversion when a reusable set resource is wanted before an
analysis. It converts standard GMT gene sets into the SNP-set block format
required by `--fastBAT-set-list`. This
is a scientific mapping step rather than a syntax-only conversion: each GMT
gene is resolved through a build-matched four-column gene list, the configured
gene window is applied, and overlapping variants are selected from a sorted
PLINK BIM file. The output copies BIM column 2 exactly, whether it contains an
rsID or another unique variant ID; it never synthesizes or translates IDs.

Chromosome normalization and the policies for duplicate GMT genes, unmapped
genes, and pathways with no BIM variants are all explicit required arguments.
The BIM must contain unique IDs and contiguous chromosome blocks with
nondecreasing positions. Variants covered by multiple genes are written once
per pathway, while `gene_variant_mapping.tsv` preserves every gene-to-variant
assignment. The utility also produces pathway summaries, unmapped-gene records,
a YAML manifest, and checksums, and publishes the directory only after every
validation succeeds.

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
  --empty-pathway-policy omit
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
- `set_annotation.oversized_set_policy`, which controls sets above GCTA's
  20,000-variant hard limit;
- `set_annotation.conversion`, including chromosome, duplicate-gene,
  unmapped-gene, empty-pathway, and output-name policies;
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
Scientific warnings identify unresolved GWAS variants and custom sets omitted
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
  --vcf study.vcf.gz \
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
  --vcf study.vcf.gz \
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
  --vcf study.vcf.gz \
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
  --vcf study.vcf.gz \
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

## Scientific sources

- [Official GCTA fastBAT documentation](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT)
- [Official GCTA mBAT-combo documentation](https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo)
- [Official GCTA fastBAT source (`sbat.cpp`)](https://github.com/jianyangqt/gcta/blob/master/main/sbat.cpp),
  which resolves summary and set entries through the reference SNP-name map;
- Bakshi et al. 2016, *Scientific Reports* 6:32894,
  [doi:10.1038/srep32894](https://doi.org/10.1038/srep32894)
- Li et al. 2023, *American Journal of Human Genetics* 110:30–43,
  [doi:10.1016/j.ajhg.2022.12.006](https://doi.org/10.1016/j.ajhg.2022.12.006)
- Benjamini and Hochberg 1995, *Journal of the Royal Statistical Society B*
  57:289–300,
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [R `p.adjust` reference implementation and definitions](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
