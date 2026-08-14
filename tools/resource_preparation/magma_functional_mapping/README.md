# MAGMA functional-mapping resources

This preparer installs pinned GRCh37/EUR inputs for conventional MAGMA,
eMAGMA, H-MAGMA, nMAGMA, and chromMAGMA. It verifies every downloaded archive,
publishes the resource tree atomically, records per-file SHA-256 checksums, and
generates a schema-validated PostGWAS configuration usable in direct and
pipeline modes.

The validated bundle contains 88 mapping definitions: 1 positional, 47 eMAGMA,
35 H-MAGMA, 4 nMAGMA, and 1 chromMAGMA.

## Prepare the resources

```bash
postgwas resources prepare magma
```

The canonical default destination is stored in the packaged resource YAML.
Override it without editing package files:

```bash
postgwas resources prepare magma \
  --output-directory /path/to/magma/functional_mapping \
  --cache-directory /path/to/download-cache
```

An existing installation is never replaced. A repeat invocation verifies the
manifest and generated PostGWAS configuration, then exits successfully.
`Rscript` is required during initial preparation to read the pinned GENCODE 26
identifier table used for eMAGMA network identifier harmonisation.

After updating only the preparer or generated-configuration schema, refresh the
metadata without downloading or modifying reference resources:

```bash
postgwas resources refresh magma
```

The command is part of the installed package, so it discovers its canonical
YAML and R helper without a repository checkout or an absolute Python-script
path. The small `prepare.py` file in this developer directory only forwards to
the same installed command.

## Use the generated configuration

```bash
postgwas magma \
  --run-config ~/Documents/software_resources/resourses/postgwas/magma/functional_mapping/configs/magma_functional_mapping.yaml \
  --magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac \
  --primary-magma-mapping positional \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --dataset-id STUDY \
  --output-directory results
```

The generated configuration selects conventional positional MAGMA by default.
Functional mappings are opt-in because each tissue/cell annotation represents a
different biological question and a separate multiple-testing family.
The same file can be supplied to `postgwas pipeline --run-config` because its
MAGMA values are stored under the canonical `modules.magma` namespace.

## Mapping contracts

- eMAGMA uses GTEx v8 tissue-specific eQTL SNP–gene assignments. PostGWAS uses
  each annotation directly and converts the optional two-column co-expression
  membership file to native MAGMA set input. The upstream annotation and network
  files use Entrez and Ensembl IDs respectively. The preparer therefore creates
  an explicit, unambiguous Ensembl–HGNC–Entrez crosswalk, preserves the original
  networks, and records per-tissue conversion/overlap counts. The configured
  positive annotation-overlap floor is an identifier-compatibility sentinel;
  it is deliberately not a high coverage requirement because significant eQTL
  genes and co-expression-network genes are tissue-specific subsets. See the
  [primary publication](https://doi.org/10.1093/bioinformatics/btab115) and
  [upstream tutorial](https://github.com/eskederks/eMAGMA-tutorial).
- H-MAGMA uses tissue/cell-specific chromatin interactions to connect noncoding
  variants to genes. See the
  [primary publication](https://doi.org/10.1038/s41593-020-0603-0) and
  [upstream repository](https://github.com/thewonlab/H-MAGMA).
- nMAGMA unions the current study's positional annotation with configured Hi-C,
  eQTL, and TOM co-expression components, removes repeated SNP assignments within
  a gene, and uses the configured `protein_coding_gene.loc` coordinates as the
  authoritative merged coordinates. Its packaged positional component uses the
  upstream zero-kilobase window. Coordinate differences in component records and
  invalid coordinate placeholders are replaced by the canonical coordinates;
  those records and genes outside the protein-coding reference are counted in
  the canonical log.
  This follows upstream `Main_script.txt` and `Merge_annot_files.R`. See the
  [primary publication](https://doi.org/10.1093/bib/bbaa298) and
  [upstream repository](https://github.com/sldrcyang/nMAGMA).
- chromMAGMA first tests configured regulatory elements with MAGMA, links them
  to genes using the supplied GeneHancer mapping, and applies the published
  lowest-element-p rule per gene. These mapped values are rankings rather than
  calibrated gene p-values, so PostGWAS does not apply ordinary gene-level
  multiple-testing correction or claim gene significance from them. See the
  [primary publication](https://doi.org/10.26508/lsa.202201446) and
  [upstream repository](https://github.com/lawrenson-lab/chromMAGMA-public).

All installed functional annotations and regulatory-element locations use
GRCh37 coordinates. External annotations use rsIDs and the generated
configuration therefore points to the official MAGMA `g1000_eur` LD reference.
PostGWAS checks exact annotation-variant overlap with BIM field 2 before an
external annotation can run.
The generated definitions also preserve the target-gene identifiers found in
the pinned inputs: Entrez for eMAGMA, Ensembl for H-MAGMA, symbols for nMAGMA,
and mixed identifiers for chromMAGMA. No method's targets are silently
relabelled.

Some pinned GitHub repositories do not contain a machine-readable licence.
Their resources are installed for local analysis with an explicit warning
in `resource_manifest.yaml`; redistribution terms must be checked separately.
