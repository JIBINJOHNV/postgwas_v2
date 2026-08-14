# Reference Resources

**All reference and resource files must be downloaded.** PostGWAS ships none of
them. Analyses require method-specific references, and compatibility is part of
the analysis—not an installation detail.

**Download link:** _to be added._

See [Resource Setup](../getting-started/resource-setup.md) for the directory
structure to create in the meantime.

## Record for every resource

Keep the resource name, upstream provider, release or version, download date,
genome build, population, chromosome convention, preparation command, checksums,
and any filters. A directory name such as `EUR` is not sufficient provenance.

## Common resource classes

- Harmonisation: reference FASTA, dbSNP or variant-identification resources,
  chromosome/build metadata, and population allele-frequency resources.
- LD annotation: population- and build-specific LD-detect region files.
- LD clumping and fine-mapping: genotype or LD references with matching variant
  identifiers, alleles, ancestry, and build.
- Imputation: the method's prepared LD and variant resources.
- LDSC: regression-weight LD scores, reference LD scores, and normally the
  HapMap3 merge-alleles list expected by the chosen LDSC resources.
- MAGMA and GCTA gene analysis: gene-location annotations plus compatible LD
  reference genotypes; gene-set analyses additionally need declared gene sets.
- PoPS and FLAMES: upstream module outputs and feature/model resources.
- MiXeR: MiXeR software resources and the LD/reference inputs expected by that
  tool.

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
