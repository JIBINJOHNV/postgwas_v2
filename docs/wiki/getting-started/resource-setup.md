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
| Positional MAGMA | Supply a matching PLINK LD reference and gene-location file | [MAGMA input requirements](../modules/magma.md#input-requirements) |
| Supported MAGMA functional mappings | Use the pinned resource downloader, then review its generated mapping configuration | [Prepare MAGMA resources](#prepare-magma-functional-mapping-resources) |
| Other downstream analyses | Supply the selected method's reference files through its public options or supported run configuration | [Reference resources by analysis](../reference/reference-resources.md#common-resource-classes) |

You may organise references under one versioned, read-only root. Harmonisation
uses `--resource-directory` or `resources.root`; other modules may take explicit
files or prefixes and do not universally expose a resource-directory option.
Record each file's source, build, population, release, and checksum.

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

## Module resources

Keep separate, clearly versioned subtrees for LD-detect BEDs, clumping LD
tables, PLINK genotype references, LDSC scores/weights, gene annotations,
MAGMA/GCTA gene sets, PoPS features, fine-mapping tools/references, MiXeR/GSA
resources, and optional FLAMES annotations. Do not reuse a resource merely
because its filename contains the expected population.

## Preflight checklist

1. Export the installed module configuration.
2. Resolve every configured path in a clean environment or container.
3. Verify build, population, chromosome names, alleles, and variant IDs.
4. Verify all chromosome files and indexes are present and non-empty.
5. Record upstream URL/version/license and checksums.
6. Test a small representative dataset and inspect join/retention rates.
7. Mount the resource tree read-only for production runs.

After these checks, follow the [Quick Start](quick-start.md). Keep the resource
manifest with the [archived analysis record](../reference/output-structure.md#archiving-a-run).
