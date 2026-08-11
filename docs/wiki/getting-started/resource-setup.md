# Resource Setup

PostGWAS does not download one universal reference bundle. Create a versioned,
read-only resource root and record where every file came from, its build and
population, preparation command, and checksum.

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

The packaged build transition is GRCh37 ↔ GRCh38. The indexed VCF comparison
sources exposed by the CLI are `1000G` and `ALFA`; ALFA EUR is the packaged
default for VCF population-frequency annotation and post-merge QC. The
separate tabular strand/MAF–EAF reference supports `ALFA`, `wgs_ukb`, `panukb`,
`1000G`, and `fingen`, and defaults to ALFA EUR. Select one source per run under
`modules.harmonisation.default_eaf`; the source name must match the filename
exactly. All five panels provide GRCh37 and GRCh38 autosomes; sex-chromosome
coverage differs by panel and build and is checked against the chromosomes in
each dataset during preflight. dbSNP source is a configuration value, not a
filename to guess.

Default AF tables use the `default_eaf` source/population and configured CHROM,
POS, ALT-as-effect, REF-as-other mapping. Comparison VCFs use `comparison_af`
and expose population INFO fields. Build-check tables use the separately
configured column mapping. Validate headers against the exported YAML before
processing a full dataset.

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
