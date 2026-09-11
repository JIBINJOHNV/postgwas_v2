# Resource Setup

Installing PostGWAS and its external programs does not prepare every analysis
reference. Select your analyses first, then acquire only the matching resources.
Keep software executables, reference datasets, and outputs in separate locations.

There is currently no universal PostGWAS reference bundle with a complete download
URL and checksum manifest. Do not assume an empty directory with the expected
names is a usable resource installation. The
[reference checklist](../reference/reference-resources.md) links each resource
class to its consuming module; the supported MAGMA downloader is described below.
When a guide does not yet provide an acquisition recipe for the required
release, that resource remains a user-supplied prerequisite, not something the
analysis command will fetch automatically.

## Choose the setup route

| Analysis route | How resources are supplied | Details |
|---|---|---|
| Raw summary statistics | Populate the configured harmonisation resource tree, including both builds used for conversion | [Harmonisation resources](#harmonisation-resource-tree) |
| Positional MAGMA | Download a matching PLINK LD reference and gene-location file; inspect the extracted companions and identifiers | [Prepare the positional example](#prepare-the-positional-magma-example) |
| Supported MAGMA functional mappings | Use the pinned resource downloader, then review its generated mapping configuration | [Prepare MAGMA resources](#prepare-magma-functional-mapping-resources) |
| Other downstream analyses | Supply the selected method's reference files through its public options or supported run configuration | [Reference resources by analysis](../reference/reference-resources.md#common-resource-classes) |

You may organise references under one versioned, read-only root. Harmonisation
uses `--resource-directory` or `resources.root`; other modules may take explicit
files or prefixes and do not universally expose a resource-directory option.
Record each file's source, build, population, release, and checksum.

## Acquire harmonisation source data

Export the installed schema-backed defaults before preparing files:

```console
postgwas config export \
  --module harmonisation \
  --style full \
  --output harmonisation.yaml
```

Review `modules.harmonisation.resource_layout`, the resource column mappings,
the selected population and build coverage. The following sources supply
ingredients, not a ready-to-use PostGWAS tree:

| Required ingredient | Authoritative acquisition route | Preparation still required |
|---|---|---|
| Assembly sequence and gene annotation | [Ensembl's GRCh37 download guide](https://grch37.ensembl.org/info/data/ftp/index.html) links FASTA and GFF3 releases; select the corresponding assembly-specific release for each build | Record the assembly and annotation release; prepare the configured chromosome FASTAs and `.fai` files and compatible GFF3. Check actual contig names, not just filenames. |
| dbSNP variants | [NCBI dbSNP VCF distribution](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/) and its release metadata | Select the correct assembly, preserve the source/version/checksums, and prepare the required chromosome files, identifier conventions and indexes. A moving `latest_release` URL is not a pinned analysis identity. |
| ALFA frequencies | [NCBI ALFA data access](https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/#ftp-download) and its release-specific files | NCBI reports these variants on GRCh38. Producing GRCh37 resources, mapping population fields and preparing PostGWAS tabular/VCF inputs are separate, validated preparation steps; renaming the source file is not build conversion. |
| Directional liftover chains | UCSC's [hg19 source directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) contains `hg19ToHg38.over.chain.gz`; its [hg38 source directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) contains `hg38ToHg19.over.chain.gz` | Choose the correct direction, review the distribution's usage terms, decompress to the configured chain path, and verify compatibility with the chosen FASTA assemblies and contig labels. Do not reverse a chain by renaming it. |
| Genome-build check tables | Build-specific variant evidence with the configured coordinate/allele columns | The repository does not currently expose a public command that downloads and assembles the complete two-build check-table bundle. Obtain or prepare documented, representative matching tables and validate their schema and coverage. |

The ALFA download is not automatically a table with PostGWAS `EUR` and other
configured column names. Preserve the population definitions and allele
orientation while preparing those fields. An external INFO reference is a
quality proxy, not a substitute for study-measured INFO.

### Existing preparation workflows and their limits

The installed, supported downloader is
[`postgwas resources prepare magma`](#prepare-magma-functional-mapping-resources).
The repository also contains specialist workflows for
[1000 Genomes Phase 3 population genotype VCFs](https://github.com/JIBINJOHNV/postgwas_v2/blob/main/tools/resource_preparation/thousand_genomes_phase3/README.md)
and [dbSNP Build 157 GRCh37 resources](https://github.com/JIBINJOHNV/postgwas_v2/blob/main/tools/resource_preparation/dbsnp_build157_grch37/README.md).
Read each workflow's README and configuration before use. Their checked-in
configurations contain machine-specific output paths that must be replaced in
a user-owned configuration, and compute/tool settings must be reviewed. The
dbSNP workflow enforces Build 157 even though its configured download URLs use
`latest_release`; verify that upstream files still match that contract.

These workflows are not a full harmonisation installer: the population-VCF
workflow does not produce every configured AF table or a ready PLINK prefix,
and the dbSNP workflow does not fill the entire per-chromosome, two-build tree.
Assembly/index preparation, frequency conversion, build-check tables, and
consistent layout still require a separately documented preparation process.
Do not start the raw-data walkthrough until that work is complete. No
end-to-end reference-preparation validation or resource-size benchmark is
claimed by this guide.

## Harmonisation resource tree

The current canonical harmonisation YAML expands these paths below
`resources.root`:

```text
<build>/default_af/tab_files/<build>_<source>_freq_chr<chromosome>.tsv.gz
<build>/external_af/vcf_files/<build>_<source>_freq_chr<chromosome>.vcf.gz
<build>/dbSNP/vcf_files/<build>_<source>_chr<chromosome>.vcf.gz
<build>/fasta_files/<build>_chr<chromosome>.fa
<build>/gff_files/<build>_ensembl.gff3.gz
chain_files/<source_build>_to_<target_build>.chain
GRCh37_38_check_files/<build>_check_file.tsv
```

Indexed files must keep their companions: `.fai` beside each FASTA, and `.tbi`
or `.csi` beside each bgzipped VCF. Supply the files required for the observed
chromosomes and for both the input and target build when conversion is required.

The packaged build transition is GRCh37 ↔ GRCh38. The indexed VCF comparison
sources exposed by the CLI are `1000G` and `ALFA`; ALFA EUR is the packaged
default for VCF population-frequency annotation and post-merge QC. The
separate tabular strand/MAF–EAF reference supports `ALFA`, `wgs_ukb`, `panukb`,
`1000G`, and `fingen`, and defaults to ALFA EUR. Select one source per run under
`modules.harmonisation.default_eaf`; the source name must match the filename
exactly. The configured panel coverage differs by chromosome and build and is
checked against each dataset during preflight; confirm that the corresponding
files actually exist in your resource release. dbSNP source is a configuration
value, not a filename to guess.

Default AF tables use the `default_eaf` source/population and configured CHROM,
POS, ALT-as-effect, REF-as-other mapping. Comparison VCFs use `comparison_af`
and expose population INFO fields. Build-check tables use the separately
configured column mapping. Validate headers against the exported YAML before
processing a full dataset.

## Prepare the positional MAGMA example

For the GRCh37/EUR example only, use the official
[MAGMA download page](https://cncr.nl/research/magma/):

1. Under **Reference data**, download the European 1000 Genomes Phase 3 archive.
   The distributed reference is Build 37; do not use it as a GRCh38 panel.
2. Under **Auxiliary files**, download the Build 37 gene-location archive. Its
   protein-coding genes use Entrez identifiers.
3. Extract each archive into a new, versioned reference directory. Locate the
   `.bed`, `.bim` and `.fam` sharing one prefix, and the `.gene.loc` file. Inspect
   the archive contents rather than assuming a tutorial placeholder is the
   downloaded filename.
4. Replace `reference/GRCh37_EUR_reference` and
   `reference/NCBI37.3.gene.loc` in the [Quick Start](quick-start.md) with those
   actual paths. Do not include `.bed` in the PLINK prefix argument.
5. Record the download source, release, checksums and applicable terms. Confirm
   build, ancestry, gene IDs and summary-statistic/BIM identifier overlap.

The `--resolve-variants-to-reference` option in the walkthrough checks exact
identifier matches; it cannot repair incompatible builds or allele conventions.
Use the [MAGMA input requirements](../modules/magma.md#input-requirements) when
selecting another ancestry, build, gene annotation or functional mapping.

## Prepare MAGMA functional-mapping resources

The public resource command currently supports the pinned MAGMA bundle. It
downloads and checksum-verifies the configured files, prepares mapping metadata,
and writes a configuration fragment for the available positional, eMAGMA,
H-MAGMA, nMAGMA, and chromMAGMA contexts:

```console
postgwas resources prepare magma \
  --output-directory reference/magma/functional_mapping
```

Inspect `postgwas resources prepare magma --help` before downloading. This
command is not a general downloader for all PostGWAS modules. Review the
generated configuration and select biologically appropriate mappings; a downloaded
tissue or cell context is not automatically appropriate for every study.
See [MAGMA](../modules/magma.md) for the distinctions between gene and regulatory
element mappings and their downstream compatibility.

To revalidate that installation later and regenerate its metadata:

```console
postgwas resources refresh magma \
  --output-directory reference/magma/functional_mapping
```

This refresh checks the installed bundle against its recorded checksums. It is
not a readiness check for your harmonisation tree or every reference used by a
downstream pipeline.

## Module resources

Keep separate, clearly versioned subtrees for LD-detect BEDs, clumping LD
tables, PLINK genotype references, LDSC scores/weights, gene annotations,
MAGMA/GCTA gene sets, PoPS features, fine-mapping tools/references, MiXeR/GSA
resources, and optional FLAMES annotations. Do not reuse a resource merely
because its filename contains the expected population.

## Preflight checklist

1. Export the installed module configuration where the command supports run YAML.
2. Resolve every configured path in a clean environment or container.
3. Verify build, population, chromosome names, alleles, and variant IDs.
4. Verify all chromosome files and indexes are present and non-empty.
5. Record upstream URL/version/license and checksums.
6. Test a small representative dataset and inspect join/retention rates.
7. Mount the resource tree read-only for production runs.

`postgwas config validate` checks configuration structure; it does not certify
every reference file. The selected analysis performs its own file and
compatibility preflight. Keep the resulting validation evidence and stop on a
missing or incompatible resource instead of treating successful CLI help or
software verification as analysis readiness.

After these checks, follow the [Quick Start](quick-start.md). Keep the resource
manifest with the [archived analysis record](../reference/output-structure.md#archiving-a-run).
