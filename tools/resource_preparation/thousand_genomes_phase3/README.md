# 1000 Genomes Phase 3 population VCF preparation

This workflow downloads every genotype-bearing VCF and its published Tabix
index from the 1000 Genomes Phase 3 `20130502` release, partitions genotypes by
the five official super-populations (`AFR`, `AMR`, `EAS`, `EUR`, and `SAS`),
recalculates genotype-derived `AC`, `AN`, `AF`, and `NS`, and writes indexed,
compressed GRCh37 VCFs.

The autosomes are concatenated into one VCF per super-population. Chromosomes
X, Y, and MT remain separate. This is required for VCF validity: the release Y
call set contains male samples only, while the MT call set is missing NA19159.
A single VCF cannot change its sample columns between contigs. No variants,
alleles, genotypes, or FORMAT fields are filtered or normalized.

The whole-genome sites-only VCF in the release directory is not downloaded. It
contains no sample genotypes and therefore cannot be partitioned by population;
its sites duplicate those represented by the genotype call sets.

## Run

The checked-in `config.yaml` is the sole source of release filenames,
population policy, paths, persistent and temporary output names, compute
settings, and cleanup policy.
It is validated strictly by Pydantic before any download or analysis begins.

Chromosome splitting and final concatenation use bounded parallel workers set
at the top of `prepare.py`. On the current 12-logical-CPU host, three concurrent
jobs combine with four `bcftools` threads per process. Restarted runs validate
and reuse each completed population/contig output independently.

```bash
python tools/resource_preparation/thousand_genomes_phase3/prepare.py \
  --config tools/resource_preparation/thousand_genomes_phase3/config.yaml \
  --stage all
```

Stages can also be run separately as `download`, `prepare`, and `validate`.
Downloads resume from `.partial` files. Existing split and final files are
reused only after their indexes and expected sample counts validate.

The release's original Tabix indexes do not contain the record-count metadata
expected by recent `bcftools index --stats`. Source validation therefore checks
the contigs recorded by the original index with `tabix --list-chroms` and
confirms the sample header. Newly generated population VCFs use current CSI
indexes and retain strict indexed-record-count validation.

Outputs include the VCFs and CSI indexes, exact group sample lists, a structured
JSON-lines run log, and a YAML resource manifest with resolved configuration,
software version, file sizes, and SHA-256 checksums.

## Scientific sources

- [Official Phase 3 release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
- [Official Phase 3 call-set README](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220)
- [Official sample/population panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel)
- [Official known-issues record](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_known_issues_20200731)

Phase 3 is GRCh37 and contains 2,504 samples from 26 populations. The panel's
`super_pop` assignments are used directly; population membership is never
inferred from sample identifiers.
