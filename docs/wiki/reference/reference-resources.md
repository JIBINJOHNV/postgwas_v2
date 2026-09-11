# Reference Resources

Analyses require method-specific reference data in addition to installed
software. The complete software installer does not supply a universal analysis
reference bundle. Some packaged methods include code or models, and the resource
command can download the supported MAGMA functional-mapping bundle; neither
replaces the other study-matched inputs required by a selected analysis.

Use [Resource Setup](../getting-started/resource-setup.md) for the harmonisation
directory layout and supported MAGMA downloader. A complete all-module reference
download URL and checksum manifest is not currently published. Where acquisition
instructions are incomplete, obtain and document the required upstream release
before running the analysis; PostGWAS does not silently substitute a panel.

## Record for every resource

Keep the resource name, upstream provider, release or version, download date,
genome build, population, chromosome convention, preparation command, checksums,
and any filters. A directory name such as `EUR` is not sufficient provenance.

## Common resource classes

| Consumer | External data or access to prepare | Detailed requirements |
|---|---|---|
| Harmonisation | Build-check tables, FASTA/index, dbSNP and frequency resources/indexes, annotation, and build-conversion chains | [Harmonisation overview](../harmonisation/overview.md), [resource tree](../getting-started/resource-setup.md#harmonisation-resource-tree) |
| LD annotation | Population- and build-specific LD-block BED files | [LD annotation](../modules/ld-annotation.md#input-requirements) |
| Standard LD clumping | Prepared manifest, forward/reverse indexed pairwise LD, and allele-aware variant inventories; not just a PLINK prefix | [LD clumping](../modules/ld-clumping.md#input-requirements) |
| Region pruning | LD-block annotations already present in the direct input VCF, or BED files for the pipeline annotation step | [LD clumping](../modules/ld-clumping.md#input-requirements) |
| COJO clumping / GCTA-COJO | Build- and ancestry-matched PLINK reference; method-specific conditioning inputs when requested | [LD clumping](../modules/ld-clumping.md), [GCTA-COJO](../../modules/gcta_cojo/README.md) |
| Fine-mapping | Matching PLINK genotypes; pipeline mode also needs the standard-clumping pairwise-LD bundle, while direct mode takes a locus file and engine-specific summary-statistic table | [Fine-mapping](../modules/fine-mapping.md#input-requirements) |
| Imputation | Prepared PRED-LD LD/variant resources plus the full harmonisation reference tree for re-harmonising imputed output | [Imputation](../modules/imputation.md#input-requirements) |
| LDSC heritability | Reference LD scores, regression-weight LD scores, and the required HapMap3 SNP/allele merge list | [LDSC heritability](../modules/ldsc.md#input-requirements) |
| MAGMA | Matching PLINK genotypes and gene mapping; matching-ID gene sets for competitive gene-set analysis | [MAGMA](../modules/magma.md#input-requirements) |
| GCTA gene/segment/set tests | Matching PLINK genotypes; gene coordinates for gene tests/GMT conversion; a native SNP-set list or GMT for set tests; fixed segments need no gene list | [GCTA gene analysis](../modules/gcta-gene.md#input-requirements) |
| MAGMAcovar / single-cell integration | Gene-property matrices, an H5AD atlas and ID mapping for scDRS, or cell-type LD scores and `.ldcts` inputs for LDSC cell-type analysis | [MAGMAcovar](../modules/magmacovar.md#input-requirements), [Single-cell](../modules/single-cell.md#input-requirements) |
| PoPS / K-POPS | Matching gene annotation and feature matrix or kernel resources | [PoPS](../modules/pops.md#input-requirements), [K-POPS](../../modules/kpops.md#inputs) |
| CALDERA / FLAMES | Compatible upstream results and method-specific annotations; selected FLAMES local modes also need VEP/CADD data | [CALDERA](../../modules/caldera.md), [FLAMES](../modules/flames.md#input-requirements) |
| MiXeR / GSA-MiXeR | Complete BIM/LD patterns; GSA annotation and GO tables for GSA analysis | [MiXeR](../modules/mixer.md#input-requirements) |
| Pathway enrichment | Gene symbols, network access, BioGRID key, registered DAVID email, and DSigDB GMT when that provider is needed | [Pathway enrichment](../modules/pathway-enrichment.md#input-requirements) |

This table is a navigation checklist, not a universal directory schema. The
selected module's guide and resolved configuration define the required file
formats, indexes, coordinate conventions, and releases. Pipeline orchestration
can produce upstream analysis artifacts; it does not download these external
references or obtain credentials.

## Compatibility checklist

Before analysis, verify:

- study and resource genome builds match;
- chromosome names and included chromosomes agree;
- effect/reference allele conventions are understood;
- the reference population is scientifically appropriate;
- required files, indexes, and chromosomes are complete;
- variant identifiers join at the expected rate;
- external-tool versions accept the prepared file formats.

Stop when these checks fail. A low join rate or allele mismatch can bias results
even when the external program exits successfully.
